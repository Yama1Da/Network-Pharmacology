"""背景互作组:状态与一键重建(Data_Management 的另一半)。"""

from fastapi import APIRouter
from fastapi import FastAPI, HTTPException, Request
from WY import config
from WY.core import jobs, transfer
from WY.core.interactome import (build_interactome, clean_old_version_caches,
                                 interactome_age_days, is_built, is_stale,
                                 require_meta)


router = APIRouter()


@router.get("/api/interactome")
def interactome_status() -> dict:
    if not is_built("Human"):
        return {"built": False, "version": config.STRING_VERSION}
    meta = require_meta("Human")
    age = interactome_age_days("Human")
    return {"built": True,
            "nodes": meta["n_nodes"], "edges": meta["n_edges"],
            "version": meta["string_version"], "confidence": meta["confidence"],
            "built_at": str(meta.get("built_at", ""))[:10],
            "age": age if age is not None else None,
            "stale": is_stale("Human")}


@router.post("/api/interactome/build")
def interactome_build() -> dict:
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

    return {"job_id": jobs.submit(task, "interactome build")}


@router.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> dict:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="no such job")
    return job
