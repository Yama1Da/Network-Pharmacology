"""来源刷新(药物/疾病管理页):单文件与批量,经 job 提交-轮询。"""

from fastapi import APIRouter
from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel
from WY.core import jobs, transfer
from WY.core.management import (enrichment_display_name, save_api_env,
                                save_ui_settings, update_all_cache_files)
from WY.core.sources import disease_spec, drug_spec


router = APIRouter()


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
    return {"job_id": job_id}
