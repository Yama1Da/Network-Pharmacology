"""来源刷新(药物/疾病管理页):单文件与批量,经 job 提交-轮询。

两类来源的刷新是相互独立的 job(独立 spec/缓存目录,update_all_cache_files
并发安全),可同时在后端各跑一个;每类在进程内记最近一次 job_id,供
/api/sources/{kind}/job 状态端点支撑页面重进恢复(同 ppi/ml/synergy 模式)。
"""

import threading
import time

from fastapi import APIRouter
from fastapi import HTTPException
from pydantic import BaseModel
from WY.core import jobs
from WY.core.management import update_all_cache_files
from WY.core.sources import disease_spec, drug_spec


router = APIRouter()

_last_refresh_job: dict[str, str | None] = {"drug": None, "disease": None}
_last_job_lock = threading.Lock()


class RefreshBody(BaseModel):
    filenames: list[str] | None = None   # None=全部更新;列表=更新选中
    hours: int | None = None             # 缓存 TTL;缺省按页面惯例(drug 12/1,disease 默认)


def _refresh_hours(kind: str, filenames: list[str] | None) -> int | None:
    if kind == "drug":
        return 12 if filenames is None else 1
    return None


@router.post("/api/sources/{kind}/refresh")
def sources_refresh(kind: str, body: RefreshBody) -> dict:
    if kind not in ("drug", "disease"):
        raise HTTPException(status_code=404, detail=f"unknown source kind '{kind}'")
    spec = drug_spec() if kind == "drug" else disease_spec()
    hours = body.hours if body.hours is not None else _refresh_hours(kind, body.filenames)
    filenames = body.filenames

    def task(progress) -> dict:
        if filenames is None:
            results = update_all_cache_files(spec, cache_hours=hours, on_progress=progress)
        else:
            results = update_all_cache_files(spec, cache_hours=hours, filenames=filenames,
                                             on_progress=progress)
        return {"rows": results.to_dict(orient="records")}

    job_id = jobs.submit(task, f"{kind} sources refresh")
    with _last_job_lock:
        _last_refresh_job[kind] = job_id
    return {"job_id": job_id}


@router.get("/api/sources/{kind}/job")
def sources_job_status(kind: str) -> dict:
    """最近一次来源刷新的在途状态(跨页离开任务照跑,页面重进恢复):
    running 照常吐 job_id 供续接;已完成但在回填窗口内(jobs.RESULT_REPLAY_WINDOW_S)
    也吐——前端 resume 命中 done 即回填结果表并刷新清单,同当场跑完的体验;
    超窗不吐,来源清单本身已是最新数据,无需回填。"""
    if kind not in ("drug", "disease"):
        raise HTTPException(status_code=404, detail=f"unknown source kind '{kind}'")
    with _last_job_lock:
        job_id = _last_refresh_job.get(kind)
    if job_id is None:
        return {"running": False}
    job = jobs.get(job_id)
    if job is None:
        return {"running": False, "job_id": job_id, "status": "unknown"}
    if job["status"] != "running" \
            and time.time() - (job.get("finished") or job["created"]) > jobs.RESULT_REPLAY_WINDOW_S:
        return {"running": False}
    return {"running": job["status"] == "running", "job_id": job_id,
            "status": job["status"], "progress": job["progress"]}
