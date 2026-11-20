"""背景互作组:状态与一键重建(Data_Management 的另一半)。"""

import threading

from fastapi import APIRouter
from fastapi import HTTPException
from WY import config
from WY.core import jobs
from WY.core.interactome import (build_interactome, clean_old_version_caches,
                                 interactome_age_days, is_built, is_stale,
                                 require_meta)


router = APIRouter()


# 在途构建登记(2026-11 审计 M:build 端点此前无防重,双击/失败重试会让两个
# job 并发下载同一批文件互写 .part/.tmp,浪费流量且必然有一边报错)。单进程
# 单用户定位下用模块级登记足够;返回同一个 job_id 让第二次点击自然并入第一次,
# 前端按既有 job_id 轮询即可(向后兼容,不认识新字段也不影响)。
_build_job_lock = threading.Lock()
_build_job_id: str | None = None


@router.get("/api/interactome")
def interactome_status() -> dict:
    if not is_built("Human"):
        status: dict = {"built": False, "version": config.STRING_VERSION}
    else:
        meta = require_meta("Human")
        age = interactome_age_days("Human")
        status = {"built": True,
                  "nodes": meta["n_nodes"], "edges": meta["n_edges"],
                  "version": meta["string_version"], "confidence": meta["confidence"],
                  "built_at": str(meta.get("built_at", ""))[:10],
                  "age": age if age is not None else None,
                  "stale": is_stale("Human")}
    # 在途构建透出:重进页面自动接上进度条(与 pk-mirror 同语义)
    with _build_job_lock:
        job_id = _build_job_id
    if job_id is not None:
        job = jobs.get(job_id)
        if job is not None and job["status"] == "running":
            status.update(building=True, job_id=job_id, progress=job["progress"])
    return status


@router.post("/api/interactome/build")
def interactome_build() -> dict:
    global _build_job_id

    def task(progress) -> dict:
        # interactome 的回调签名是 (stage, fraction, detail),适配成 job 进度的 (fraction, text);
        # 直接透传会因参数个数不符每次抛 TypeError,前端全程看不到构建进度
        stage_text = {
            "download_links": "下载 STRING links",
            "download_info": "下载 STRING info",
            "mirror_skip": "跳过不可用镜像",
            "fallback": "并发下载未完成,改用单流(约需 10 分钟以上)",
            "processing": "解析交互组",
            "writing": "写入缓存",
        }

        def cb(stage, fraction=None, detail=""):
            text = stage_text.get(stage, stage)
            if detail == "cached":
                text = f"{text}(已有缓存)"
            elif detail:
                text = f"{text} · {detail}"
            progress(fraction=fraction, text=text)

        meta_new = build_interactome(progress_cb=cb)
        cleaned = clean_old_version_caches("Human")
        return {"nodes": meta_new["n_nodes"], "edges": meta_new["n_edges"],
                "cleaned": cleaned}

    # 查在途 + 提交 + 登记必须同锁原子完成:分两段的话两个并发请求会都看到
    # "无在途"而各自提交,防重形同虚设。submit 只登记并起线程(无 I/O),持锁
    # 代价可忽略;job 已结束(done/error)则清旧登记照常提交新构建。
    with _build_job_lock:
        if _build_job_id is not None:
            job = jobs.get(_build_job_id)
            if job is not None and job["status"] == "running":
                return {"job_id": _build_job_id, "already_running": True}
        job_id = jobs.submit(task, "interactome build")
        _build_job_id = job_id
    return {"job_id": job_id}


@router.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> dict:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="no such job")
    return job
