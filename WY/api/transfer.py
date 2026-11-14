"""数据迁移:设置页的结果导出/两段式导入。"""

from fastapi import APIRouter
import os
import uuid
from fastapi import HTTPException, Request
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse
from pydantic import BaseModel
from WY import config
from WY.core import jobs, transfer
from .common import (_settings_snapshot)


router = APIRouter()


# 业务逻辑(打包范围/manifest 校验/解压安检/冲突落盘)都在 core.transfer,
# 这里只做 HTTP 形态:job 提交、文件流式收发、会话生命周期、错误码转换。

class TransferExportBody(BaseModel):
    scopes: list[str]


@router.post("/api/transfer/export")
def transfer_export(body: TransferExportBody) -> dict:
    """导出为 zip:走 job 打包(core.build_export_bundle),完成后经下载端点取包。"""
    unknown = sorted({s for s in body.scopes if s not in transfer.ALL_SCOPES})
    if unknown:
        raise HTTPException(status_code=400, detail=f"unknown scopes: {unknown}")
    scopes = [s for s in transfer.ALL_SCOPES if s in body.scopes]
    if not scopes:
        raise HTTPException(status_code=400, detail="empty scopes")
    # 空库快速失败:免得 job 打包出一个只有 manifest 的空包
    if not any(s == "ui_settings" or transfer.list_scope_items(s) for s in scopes):
        raise HTTPException(status_code=400,
                            detail="nothing to export in selected scopes")

    def task(progress) -> dict:
        transfer.cleanup_transfer_staging()
        os.makedirs(config.TRANSFER_DIR, exist_ok=True)
        out = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
        return transfer.build_export_bundle(scopes, out, progress=progress)

    job_id = jobs.submit(task, "results export")
    return {"job_id": job_id}


@router.get("/api/transfer/export/{filename}")
def transfer_export_download(filename: str) -> FileResponse:
    """取回导出包:文件名白名单(basename + 固定命名形态)防路径探测。"""
    if not transfer.EXPORT_FILE_RE.fullmatch(filename):
        raise HTTPException(status_code=400, detail=f"invalid export filename: {filename!r}")
    path = os.path.join(config.TRANSFER_DIR, filename)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="export bundle not found (expired?)")
    return FileResponse(
        path, media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'})


@router.post("/api/transfer/import")
async def transfer_import(request: Request) -> dict:
    """两段式导入·第一步:zip 原始请求体流式落盘(不引 multipart 依赖)→
    解压暂存 + 安检 + 冲突分析,返回 (sid, 预览报告) 供前端逐条决策。"""
    transfer.cleanup_transfer_staging()
    os.makedirs(config.TRANSFER_DIR, exist_ok=True)
    upload = os.path.join(config.TRANSFER_DIR, f"upload_{uuid.uuid4().hex[:12]}.zip")
    received = 0
    try:
        with open(upload, "wb") as f:
            async for chunk in request.stream():
                received += len(chunk)
                if received > config.TRANSFER_MAX_UPLOAD_BYTES:
                    raise HTTPException(status_code=413,
                                        detail="upload exceeds transfer size limit")
                f.write(chunk)
    except BaseException:
        transfer.remove_quietly(upload)
        raise
    if received == 0:
        transfer.remove_quietly(upload)
        raise HTTPException(status_code=400, detail="empty upload")
    try:
        # 解压安检放线程池:大包秒级耗时,同步做会阻塞事件循环拖慢其他请求
        sid, report = await run_in_threadpool(transfer.create_import_session, upload)
    except transfer.TransferError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return {"sid": sid, "report": report}


class TransferDecision(BaseModel):
    kind: str
    name: str
    action: str            # overwrite | skip | rename(ui_settings 只认前两种的 overwrite/skip)
    new_name: str = ""     # action=rename 时的用户自填新名


class TransferCommitBody(BaseModel):
    decisions: list[TransferDecision]


@router.post("/api/transfer/import/{sid}/commit")
def transfer_import_commit(sid: str, body: TransferCommitBody) -> dict:
    """两段式导入·第二步:按用户决策落盘。存在失败项时保留会话以便调整重试,
    全部落定(含跳过)即销毁会话与暂存。"""
    session = transfer.get_import_session(sid)
    if session is None:
        raise HTTPException(status_code=404,
                            detail="import session expired, upload the bundle again")
    try:
        report = transfer.commit_import(
            session["staging"], [d.model_dump() for d in body.decisions])
    except transfer.TransferError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    ui_snapshot = None
    if any(a["kind"] == "ui_settings" for a in report["applied"]):
        ui_snapshot = _settings_snapshot()
    if not report["failed"]:
        transfer.drop_import_session(sid)
    return {"report": report, "ui_snapshot": ui_snapshot}


@router.delete("/api/transfer/import/{sid}")
def transfer_import_cancel(sid: str) -> dict:
    """放弃导入:销毁会话并清空暂存目录。"""
    if transfer.get_import_session(sid) is None:
        raise HTTPException(status_code=404, detail="import session not found")
    transfer.drop_import_session(sid)
    return {"cancelled": True}
