"""PK-DB 镜像:状态与一键构建/更新(数据管理页"PK-DB 镜像"卡片的后端)。

交互照 api/interactome.py 的防重模式:单进程单用户,双击/失败重试的并发提交
会互相覆盖镜像文件——"查在途 + 提交 + 登记"同锁原子完成,返回同一 job_id
让第二次点击自然并入第一次,前端按既有 job_id 轮询即可。
"""

import logging
import threading

from fastapi import APIRouter
from pydantic import BaseModel
from WY.core import jobs
from WY.core.pk_mirror import build_mirror, mirror_status
from WY.core.reachability import normalization_status

logger = logging.getLogger(__name__)

router = APIRouter()


_build_job_lock = threading.Lock()
_build_job_id: str | None = None


class BuildBody(BaseModel):
    force: bool = False   # True=已存在的研究也重拉(管理页"更新"按钮;失败清单随新构建刷新)


@router.get("/api/pk-mirror")
def pk_mirror_status() -> dict:
    """磁盘镜像状态 + 归一化缓存摘要 + 在途构建(重进页面自动接上进度条)。
    normalization 块由 reachability 层供给(镜像层不回头读下游产物,分层红线)。"""
    status = mirror_status()
    try:
        status["normalization"] = normalization_status()
    except Exception:
        logger.exception("归一化状态读取失败")
        status["normalization"] = None
    with _build_job_lock:
        job_id = _build_job_id
    if job_id is not None:
        job = jobs.get(job_id)
        if job is not None and job["status"] == "running":
            status.update(building=True, job_id=job_id, progress=job["progress"])
    return status


@router.post("/api/pk-mirror/build")
def pk_mirror_build(body: BuildBody | None = None) -> dict:
    global _build_job_id
    force = bool(body.force) if body is not None else False

    def task(progress) -> dict:
        result = build_mirror(force=force, progress=progress)
        # 级联归一化(2026-11-25 负责人定稿,方案 §3.1 分层的用意所在):
        # 镜像落盘后顺手跑三件套(force=True——数据随镜像刷新,规则版本不变),
        # 本地 5 秒级;新部署不再出现"镜像建了、可达性列全 —"的断链。
        # 失败不推翻镜像结果:显式记 normalize_error,可再次"更新"重试。
        logger.debug("PK 镜像构建完成(studies=%s),开始归一化级联",
                     result.get("studies"))
        try:
            from WY.core.reachability import build_normalized_cache
            norm_paths = build_normalized_cache(force=True)
            result["normalized"] = True
            logger.debug("归一化级联完成: %s", norm_paths)
        except Exception:
            result["normalized"] = False
            result["normalize_error"] = "normalize_failed"
            logger.exception("归一化级联失败(镜像本身已落盘;再次'更新'可重试)")
        return result

    with _build_job_lock:
        if _build_job_id is not None:
            job = jobs.get(_build_job_id)
            if job is not None and job["status"] == "running":
                return {"job_id": _build_job_id, "already_running": True}
        job_id = jobs.submit(task, "PK-DB mirror build")
        _build_job_id = job_id
    return {"job_id": job_id}
