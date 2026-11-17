"""PK-DB 镜像:状态与一键构建/更新(数据管理页"PK-DB 镜像"卡片的后端)。

交互照 api/interactome.py 的防重模式:单进程单用户,双击/失败重试的并发提交
会互相覆盖镜像文件——"查在途 + 提交 + 登记"同锁原子完成,返回同一 job_id
让第二次点击自然并入第一次,前端按既有 job_id 轮询即可。
"""

import threading

from fastapi import APIRouter
from pydantic import BaseModel
from WY.core import jobs
from WY.core.pk_mirror import build_mirror, mirror_status

router = APIRouter()


_build_job_lock = threading.Lock()
_build_job_id: str | None = None


class BuildBody(BaseModel):
    force: bool = False   # True=已存在的研究也重拉(管理页"更新"按钮;失败清单随新构建刷新)


@router.get("/api/pk-mirror")
def pk_mirror_status() -> dict:
    return mirror_status()


@router.post("/api/pk-mirror/build")
def pk_mirror_build(body: BuildBody | None = None) -> dict:
    global _build_job_id
    force = bool(body.force) if body is not None else False

    def task(progress) -> dict:
        return build_mirror(force=force, progress=progress)

    with _build_job_lock:
        if _build_job_id is not None:
            job = jobs.get(_build_job_id)
            if job is not None and job["status"] == "running":
                return {"job_id": _build_job_id, "already_running": True}
        job_id = jobs.submit(task, "PK-DB mirror build")
        _build_job_id = job_id
    return {"job_id": job_id}
