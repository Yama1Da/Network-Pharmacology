"""数据迁移:结果文件的导出打包与两段式导入(设置页「数据迁移」板块的后端)。

导出只打四类用户数据:PPI+富集等分析结果目录(downloads 根)、分子对接 star 收藏、
GEO 筛选结果、可选的 ui_settings.json。cache/日志一律不进包——密钥绝不进包
(密钥只存在于 API.env);打包时再跳过名为 API.env/*.env 的文件作双保险。

导入两段式:上传的 zip 先解压到暂存目录做分析(manifest 校验/条目路径安检/
zip 炸弹防护/与现有数据的冲突检测),前端逐条让用户决策(覆盖/跳过/改名)后才落盘。
导入会话仿 docking_session:进程内 dict + TTL 清扫,暂存目录随会话过期一并删除;
进程重启后内存登记丢失,残留暂存由启动清扫按目录年龄兜底删除。
"""

import json
import os
import re
import shutil
import time
import uuid
import zipfile
from datetime import datetime

# 导出包名与包内 manifest 的固定形态;导入端按它识别"本程序导出的包"
EXPORT_FILE_RE = re.compile(r"^npap_export_\d{8}_\d{6}\.zip$")
MANIFEST_NAME = "manifest.json"
MANIFEST_VERSION = 1
# 目录集合类范围;ui_settings 是单文件范围,单独处理
DIR_KINDS = ("analysis", "star", "geo")
ALL_SCOPES = (*DIR_KINDS, "ui_settings")
# analysis 范围排除 downloads 根下的三个产物区目录(各自是独立范围);
# 目录名与 config.py 中 MD_RESULT_DIR/MD_STAR_DIR/GEO_RESULT_DIR 的字面量一致,写死防漂移
_ANALYSIS_EXCLUDES = ("MDresult", "star", "GEOresult")

_SESSIONS: dict[str, dict] = {}


class TransferError(ValueError):
    """导入包不合法(校验不通过);消息直接面向前端展示,api 层转 400。"""


def remove_quietly(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass


# ======== 范围与路径解析 ========

def _scope_root(kind: str) -> str:
    """范围 kind → 磁盘根目录;导入落盘与导出收集共用同一映射。"""
    from WY import config

    if kind == "analysis":
        return config.DOWNLOADS_DIR
    if kind == "star":
        return config.MD_STAR_DIR
    if kind == "geo":
        return config.GEO_RESULT_DIR
    raise TransferError(f"unknown scope kind: {kind!r}")


def list_scope_items(kind: str) -> list[str]:
    """列出某范围内可导出的条目名(仅目录,隐藏项跳过,与 files._list_to_table 同口径)。"""
    root = _scope_root(kind)
    if not os.path.isdir(root):
        return []
    excludes = _ANALYSIS_EXCLUDES if kind == "analysis" else ()
    return sorted(
        n for n in os.listdir(root)
        if not n.startswith(".") and n not in excludes
        and os.path.isdir(os.path.join(root, n))
    )


def _ensure_item_name(name) -> str:
    """条目名必须是单段安全路径组件:防 manifest 伪造路径(../、盘符、分隔符)。"""
    if not isinstance(name, str) or not name or name != os.path.basename(name) \
            or name in (".", "..") or "\\" in name or "/" in name or ":" in name:
        raise TransferError(f"illegal item name in bundle: {name!r}")
    return name


# ======== 导出 ========

def export_filename() -> str:
    return "npap_export_" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".zip"


def _manifest() -> dict:
    """manifest 头部(不含 items):schema 版本随包走,导入端比对后决定是否提示。"""
    from WY import config

    try:
        from importlib.metadata import version
        app_version = version("npap")
    except Exception:
        app_version = "unknown"
    return {
        "app": "npap",
        "manifest_version": MANIFEST_VERSION,
        "app_version": app_version,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "schema": {
            "target_filter_version": config.TARGET_FILTER_VERSION,
            "ml_screening_version": config.ML_SCREENING_VERSION,
            "string_version": config.STRING_VERSION,
        },
    }


def _is_env_file(name: str) -> bool:
    """密钥文件双保险:打包与导入都不碰 API.env/*.env(密钥绝不进包)。"""
    base = os.path.basename(name).lower()
    return base == "api.env" or base.endswith(".env")


def build_export_bundle(scopes, out_path: str, progress=None) -> dict:
    """把选定范围打包为 zip(根放 manifest,可选 ui_settings.json,条目在 items/<kind>/下)。

    progress(fraction, text) 供 job 进度回调;返回统计供 job result 与前端提示。"""
    from WY import config

    scopes = tuple(scopes)
    unknown = [s for s in scopes if s not in ALL_SCOPES]
    if unknown:
        raise TransferError(f"unknown scopes: {unknown}")
    if not scopes:
        raise TransferError("empty scopes")

    dir_scopes = [k for k in DIR_KINDS if k in scopes]
    planned = {kind: list_scope_items(kind) for kind in dir_scopes}
    has_ui = "ui_settings" in scopes and os.path.isfile(config.UI_SETTINGS_FILE)
    if not any(planned.values()) and not has_ui:
        raise TransferError("nothing to export in selected scopes")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    items: list[dict] = []
    ui_exported = False
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:

        def add_file(path: str, arcname: str) -> None:
            if _is_env_file(arcname):
                return  # 双保险:密钥文件绝不进包
            zf.write(path, arcname)

        for i, kind in enumerate(dir_scopes):
            names = planned[kind]
            for name in names:
                root = os.path.join(_scope_root(kind), name)
                for r, _, files in os.walk(root):
                    for fn in files:
                        fp = os.path.join(r, fn)
                        rel = os.path.relpath(fp, root).replace(os.sep, "/")
                        add_file(fp, f"items/{kind}/{name}/{rel}")
                items.append({"kind": kind, "name": name})
            if progress is not None:
                progress(fraction=(i + 1) / (len(dir_scopes) + 1),
                         text=f"packing {kind}: {len(names)} item(s)")

        if "ui_settings" in scopes and os.path.isfile(config.UI_SETTINGS_FILE):
            add_file(config.UI_SETTINGS_FILE, "ui_settings.json")
            ui_exported = True

        manifest = _manifest()
        manifest["items"] = items
        manifest["ui_settings"] = ui_exported
        zf.writestr(MANIFEST_NAME, json.dumps(manifest, ensure_ascii=False, indent=2))

    counts = {kind: sum(1 for it in items if it["kind"] == kind) for kind in DIR_KINDS}
    return {"filename": os.path.basename(out_path), "size": os.path.getsize(out_path),
            "counts": counts, "ui_settings": ui_exported}


# ======== 导入·分析与落盘 ========

def _safe_target(root: str, arcname: str) -> tuple[str, list[str]]:
    """zip 条目名 → (root 下的落盘路径, 路径段);越界条目一律拒收。

    不用 zf.extract/extractall:它们对 ".." 是静默剥除,这里选择显式拒绝整个包,
    保证"包有问题"尽早暴露而不是被悄悄修好。"""
    if "\\" in arcname:
        raise TransferError(f"illegal entry name (backslash): {arcname!r}")
    normalized = arcname.replace("\\", "/")
    if normalized.startswith("/"):
        raise TransferError(f"absolute path in bundle: {arcname!r}")
    if os.path.splitdrive(normalized)[0]:
        raise TransferError(f"drive-qualified path in bundle: {arcname!r}")
    parts = normalized.rstrip("/").split("/")
    if any(p in ("", ".", "..") or ":" in p for p in parts):
        raise TransferError(f"illegal path component in bundle entry: {arcname!r}")
    target = os.path.join(root, *parts)
    root_abs = os.path.abspath(root)
    target_abs = os.path.abspath(target)
    if os.path.commonpath([root_abs, target_abs]) != root_abs:
        raise TransferError(f"entry escapes staging root: {arcname!r}")
    return target, parts


def _junk_entry(parts: list[str]) -> bool:
    """打包工具与操作系统塞进来的杂项(MacOS 元数据/缩略图缓存),解压时静默跳过。"""
    base = parts[-1]
    return parts[0] == "__MACOSX" or base in (".DS_Store", "Thumbs.db") \
        or base.startswith("._")


def _validate_manifest(manifest) -> list[dict]:
    """manifest 结构校验:版本不超前、条目 kind 合法、条目名是安全路径组件、无重复。"""
    if not isinstance(manifest, dict):
        raise TransferError("manifest.json is not a JSON object")
    mv = manifest.get("manifest_version")
    if not isinstance(mv, int) or isinstance(mv, bool):
        raise TransferError("manifest.json missing integer manifest_version")
    if mv > MANIFEST_VERSION:
        raise TransferError(
            f"bundle manifest v{mv} is newer than supported v{MANIFEST_VERSION}; "
            "please upgrade this app first")
    items = manifest.get("items")
    if not isinstance(items, list):
        raise TransferError("manifest.json missing items list")
    seen = set()
    for it in items:
        if not isinstance(it, dict):
            raise TransferError("manifest items must be objects")
        kind = it.get("kind")
        if kind not in DIR_KINDS:
            raise TransferError(f"unknown item kind in manifest: {kind!r}")
        name = _ensure_item_name(it.get("name"))
        if (kind, name) in seen:
            raise TransferError(f"duplicate item in manifest: {kind}/{name}")
        seen.add((kind, name))
    return manifest


def _read_manifest(zf: zipfile.ZipFile) -> dict:
    try:
        return json.loads(zf.read(MANIFEST_NAME).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TransferError("manifest.json is not valid JSON") from exc


def _copy_capped(src, dst, cap: int) -> None:
    """限长拷贝:最多写 cap+1 字节,写满后仍探得到数据即判定谎报
    (zipfile 自身对"声明大小≠实际"会先炸 CRC,这里是不依赖其行为的兜底)。"""
    remaining = cap + 1
    while True:
        n = min(256 * 1024, remaining)
        chunk = src.read(n if n > 0 else 1)  # 写满限额后仍探 1 字节,不留静默截断
        if not chunk:
            return
        remaining -= len(chunk)
        if remaining < 0:
            raise TransferError("entry larger than declared size: possible zip bomb")
        dst.write(chunk)


def _geo_item_invalid(name: str, run_dir: str) -> str | None:
    """geo 条目有效返回 None,否则返回不可导入的原因(无效条目放进去也不会出现在历史列表)。"""
    from WY.core.ml_results import run_tag

    if run_tag(name) is None:
        return "directory name lacks _YYYYMMDD_HHMMSS suffix"
    if not os.path.isfile(os.path.join(run_dir, "summary.json")):
        return "missing summary.json"
    return None


def _schema_warnings(schema: dict) -> list[str]:
    """导出包与当前程序的算法/数据版本差异:只提示不阻断(结果文件是静态产物)。"""
    from WY import config

    warnings = []
    if not isinstance(schema, dict):
        return warnings
    tfv = schema.get("target_filter_version")
    mlv = schema.get("ml_screening_version")
    sv = schema.get("string_version")
    if isinstance(tfv, int) and tfv != config.TARGET_FILTER_VERSION:
        warnings.append(
            f"bundle was built with target filter v{tfv}, this app uses "
            f"v{config.TARGET_FILTER_VERSION}; analysis results are for reference only")
    if isinstance(mlv, int) and mlv != config.ML_SCREENING_VERSION:
        warnings.append(
            f"bundle was built with ML screening v{mlv}, this app uses "
            f"v{config.ML_SCREENING_VERSION}; imported GEO runs will be marked stale")
    if isinstance(sv, str) and sv != config.STRING_VERSION:
        warnings.append(
            f"bundle was built with STRING {sv}, this app uses {config.STRING_VERSION}; "
            "static results are unaffected")
    return warnings


def analyze_bundle(zip_path: str, staging_dir: str) -> dict:
    """两段式导入·第一步:校验并解压导入包到暂存目录,返回给前端的预览报告。"""
    from WY import config

    try:
        zf = zipfile.ZipFile(zip_path)
    except zipfile.BadZipFile as exc:
        raise TransferError("not a valid zip file") from exc

    with zf:
        infos = zf.infolist()
        if MANIFEST_NAME not in {i.filename for i in infos}:
            raise TransferError("missing manifest.json: not an npap export bundle")
        manifest = _validate_manifest(_read_manifest(zf))

        # zip 炸弹防护:先按中央目录声明的未压缩总量设限,再动手解压
        total = sum(i.file_size for i in infos)
        if total > config.TRANSFER_MAX_EXTRACT_BYTES:
            raise TransferError(
                f"bundle decompressed size ({total} bytes) exceeds the limit; "
                "possible zip bomb")

        # 条目必须都能对上 manifest(或包级白名单文件);先全量安检,再动手解压
        manifest_items = {(it["kind"], it["name"]) for it in manifest["items"]}
        ui_in_bundle = manifest.get("ui_settings") is True
        plan: list[tuple[zipfile.ZipInfo, str, list[str]]] = []
        for info in infos:
            target, parts = _safe_target(staging_dir, info.filename)
            if _junk_entry(parts):
                continue
            if info.filename == MANIFEST_NAME:
                pass
            elif info.filename == "ui_settings.json":
                if not ui_in_bundle:
                    raise TransferError("unexpected ui_settings.json in bundle")
            elif len(parts) >= 4 and parts[0] == "items" \
                    and (parts[1], parts[2]) in manifest_items:
                pass
            else:
                raise TransferError(f"unexpected entry in bundle: {info.filename!r}")
            plan.append((info, target, parts))

        os.makedirs(staging_dir, exist_ok=True)
        # manifest 登记的每个条目先把目录立起来:空目录(导出时 zip 里没有任何文件)
        # 也要能原样导回,否则 commit 会因暂存缺目录而报 not found
        for kind, name in manifest_items:
            os.makedirs(os.path.join(staging_dir, "items", kind, name), exist_ok=True)
        stat: dict[tuple[str, str], dict] = {}
        try:
            for info, target, parts in plan:
                if info.is_dir():
                    os.makedirs(target, exist_ok=True)
                    continue
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with zf.open(info) as src, open(target, "wb") as dst:
                    _copy_capped(src, dst, info.file_size)
                if len(parts) >= 4 and parts[0] == "items":
                    s = stat.setdefault((parts[1], parts[2]), {"files": 0, "bytes": 0})
                    s["files"] += 1
                    s["bytes"] += info.file_size
        except TransferError:
            raise
        except Exception as exc:
            # 中途损坏/系统层失败(截断包 BadZipFile、CRC 不符、Windows 路径超长等):
            # 归一为 TransferError,api 层转 400,而不是裸异常变 500
            raise TransferError(f"bundle is corrupted or unreadable: {exc}") from exc

    return _build_report(manifest, stat, staging_dir)


def _build_report(manifest: dict, stat: dict, staging_dir: str) -> dict:
    """从 manifest 与解压统计组装预览报告:条目清单/冲突标记/geo 有效性/版本差异提示。"""
    from WY import config

    items = []
    for it in manifest["items"]:
        kind, name = it["kind"], it["name"]
        s = stat.get((kind, name), {"files": 0, "bytes": 0})
        entry = {"kind": kind, "name": name, "files": s["files"], "bytes": s["bytes"],
                 "conflict": os.path.isdir(os.path.join(_scope_root(kind), name))}
        if kind == "geo":
            reason = _geo_item_invalid(
                name, os.path.join(staging_dir, "items", kind, name))
            if reason:
                entry["invalid"] = True
                entry["invalid_reason"] = reason
        items.append(entry)

    ui = {"present": manifest.get("ui_settings") is True}
    if ui["present"]:
        path = os.path.join(staging_dir, "ui_settings.json")
        from WY.core.caching import load_json_file
        raw = load_json_file(path)
        ui["keys"] = len(_filter_ui_settings(raw)) if isinstance(raw, dict) else 0
        ui["conflict"] = os.path.isfile(config.UI_SETTINGS_FILE)
    return {
        "manifest": {"app_version": manifest.get("app_version", ""),
                     "created_at": manifest.get("created_at", ""),
                     "schema": manifest.get("schema", {})},
        "items": items,
        "ui_settings": ui,
        "warnings": _schema_warnings(manifest.get("schema", {})),
    }


def _filter_ui_settings(raw: dict) -> dict:
    """包内 ui_settings 按白名单+类型过滤(校验语义与 management.save_ui_settings /
    config.apply_ui_settings_overrides 一致),坏值与未知键直接丢弃:设置文件可以被
    整体替换,但绝不写入白名单外的内容。"""
    from WY import config

    out = {}
    for k, v in raw.items():
        if k == "CHEMBL_TARGET_BLACKLIST":
            if isinstance(v, list) and all(isinstance(i, str) for i in v):
                out[k] = v
        elif k == "STRING_VERSION":
            sv = config.normalized_string_version(v)
            if sv:
                out[k] = sv
        elif k == "UI_FONT_SCALE":
            if isinstance(v, (int, float)) and not isinstance(v, bool) \
                    and 0.8 <= v <= 1.5:
                out[k] = float(v)
        elif k == "UI_THEME":
            if isinstance(v, str) and v.strip().lower() in ("light", "dark"):
                out[k] = v.strip().lower()
        elif k == "UI_LANG":
            if isinstance(v, str) and v.strip().lower() in config.UI_LANG_VALUES:
                out[k] = v.strip().lower()
        elif k == "LOG_LEVEL":
            if isinstance(v, str) and v.strip().upper() in (
                    "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"):
                out[k] = v.strip().upper()
        elif k in config.UI_SETTINGS_KEYS:
            if isinstance(v, bool):
                out[k] = v
            elif isinstance(v, (int, float)) and not isinstance(v, bool):
                out[k] = v
    return out


def commit_import(staging_dir: str, decisions: list[dict]) -> dict:
    """两段式导入·第二步:按用户决策落盘。skip 不动,overwrite 先隔离旧数据再移入,
    rename 校验新名后移入;ui_settings 只认 覆盖/跳过。"""
    applied, skipped, failed = [], [], []
    for d in decisions:
        kind, name, action = d.get("kind"), d.get("name"), d.get("action")
        try:
            if action == "skip":
                skipped.append({"kind": str(kind), "name": str(name)})
                continue
            if kind == "ui_settings":
                if action != "overwrite":
                    raise TransferError(f"unknown action for ui_settings: {action!r}")
                applied.append(_apply_ui_settings(staging_dir))
                continue
            if kind not in DIR_KINDS:
                raise TransferError(f"unknown kind: {kind!r}")
            _ensure_item_name(name)
            src = os.path.join(staging_dir, "items", kind, name)
            if not os.path.isdir(src):
                raise TransferError(f"not found in bundle: {kind}/{name}")
            if kind == "geo":
                reason = _geo_item_invalid(name, src)
                if reason:
                    raise TransferError(f"invalid geo run ({reason})")
            if action == "overwrite":
                applied.append(_move_into_place(staging_dir, kind, name, name))
            elif action == "rename":
                new_name = _ensure_item_name((d.get("new_name") or "").strip())
                if os.path.isdir(os.path.join(_scope_root(kind), new_name)):
                    raise TransferError(f"new name already exists: {new_name}")
                applied.append(_move_into_place(staging_dir, kind, name, new_name))
            else:
                raise TransferError(f"unknown action: {action!r}")
        except (TransferError, OSError) as exc:
            failed.append({"kind": str(kind), "name": str(name), "reason": str(exc)})
    return {"applied": applied, "skipped": skipped, "failed": failed}


def _move_into_place(staging_dir: str, kind: str, name: str, dest_name: str) -> dict:
    """把暂存条目移入目标位置。目标已存在时(覆盖语义)旧目录先 rename 进暂存隔离区
    (同卷瞬时),新内容落位成功后再删;落位失败把隔离区移回,保证旧数据不丢。"""
    src = os.path.join(staging_dir, "items", kind, name)
    dest = os.path.join(_scope_root(kind), dest_name)
    quarantine = None
    if os.path.isdir(dest):
        quarantine = os.path.join(staging_dir, "_replaced", dest_name)
        os.makedirs(os.path.dirname(quarantine), exist_ok=True)
        if os.path.exists(quarantine):
            shutil.rmtree(quarantine, ignore_errors=True)
        shutil.move(dest, quarantine)
    try:
        shutil.move(src, dest)
    except OSError:
        if quarantine is not None and not os.path.isdir(dest):
            shutil.move(quarantine, dest)  # 回滚:旧数据物归原位
        raise
    if quarantine is not None:
        shutil.rmtree(quarantine, ignore_errors=True)
    return {"kind": kind, "name": dest_name}


def _apply_ui_settings(staging_dir: str) -> dict:
    """用包内 ui_settings.json 整体替换本地设置(白名单过滤),并回填运行时。"""
    from WY import config
    from WY.core.caching import load_json_file, save_json_file

    raw = load_json_file(os.path.join(staging_dir, "ui_settings.json"))
    if not isinstance(raw, dict):
        raise TransferError("ui_settings.json in bundle is not a JSON object")
    cleaned = _filter_ui_settings(raw)
    save_json_file(config.UI_SETTINGS_FILE, cleaned)
    config.apply_ui_settings_overrides()
    return {"kind": "ui_settings", "name": "ui_settings.json", "keys": sorted(cleaned)}


# ======== 导入会话(进程内 dict + TTL,暂存目录随会话清理) ========

def create_import_session(zip_path: str) -> tuple[str, dict]:
    """解压上传包并登记导入会话,返回 (sid, 预览报告)。上传的 zip 解压后即删。"""
    from WY import config

    cleanup_transfer_staging()
    sid = uuid.uuid4().hex[:12]
    staging = os.path.join(config.TRANSFER_DIR, f"import_{sid}")
    try:
        report = analyze_bundle(zip_path, staging)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)  # 失败不留半截暂存等 TTL
        raise
    finally:
        remove_quietly(zip_path)  # 无论成败,zip 都已不再需要
    _SESSIONS[sid] = {"_touched": time.time(), "staging": staging, "report": report}
    return sid, report


def get_import_session(sid: str) -> dict | None:
    session = _SESSIONS.get(sid)
    if session is not None:
        session["_touched"] = time.time()  # 摸一次续命:按最后访问计时
    return session


def drop_import_session(sid: str) -> None:
    session = _SESSIONS.pop(sid, None)
    if session:
        shutil.rmtree(session.get("staging", ""), ignore_errors=True)


def cleanup_transfer_staging() -> None:
    """清扫 TRANSFER_DIR:过期会话(内存登记)、孤儿暂存与半截上传(进程重启残留)、
    过期导出包。新建会话与 api 启动时都会调用。"""
    from WY import config

    now = time.time()
    for sid in [s for s, d in _SESSIONS.items()
                if d.get("_touched", 0) < now - config.TRANSFER_SESSION_TTL]:
        drop_import_session(sid)
    if not os.path.isdir(config.TRANSFER_DIR):
        return
    live = {os.path.basename(d.get("staging", "")) for d in _SESSIONS.values()}
    for name in os.listdir(config.TRANSFER_DIR):
        path = os.path.join(config.TRANSFER_DIR, name)
        try:
            age = now - os.path.getmtime(path)
        except OSError:
            continue
        if name.startswith("import_") and name not in live \
                and age > config.TRANSFER_SESSION_TTL:
            shutil.rmtree(path, ignore_errors=True)
        elif name.startswith("upload_") and age > config.TRANSFER_SESSION_TTL:
            remove_quietly(path)
        elif EXPORT_FILE_RE.fullmatch(name) and age > config.TRANSFER_EXPORT_TTL:
            remove_quietly(path)
