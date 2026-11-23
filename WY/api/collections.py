"""通用文件集合:六类管理页的清单/打包/删除。"""

from fastapi import APIRouter
import os
import shutil
from fastapi import HTTPException
from fastapi.responses import Response
from pydantic import BaseModel
from WY import config
from WY.core.enrichment import parse_enrichment_filename
from WY.core.files import (list_cache_files, list_download_dirs,
                           make_dirs_zip, make_ppi_zip)
from WY.core.management import enrichment_display_name
from WY.core.sources import disease_spec, drug_spec
from WY.core.sources import file_path as spec_file_path
from WY.core.sources import merged_status
from .common import (_ZIP_NAMES, _zip_response)


router = APIRouter()


# mode "files" 条目为文件(os.remove 删除,make_ppi_zip 打包);
# mode "dirs" 条目为目录(rmtree 删除,make_dirs_zip 打包)。

def _list_mtime(path: str) -> str:
    import datetime

    try:
        return datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime(
            "%Y-%m-%d %H:%M:%S")
    except OSError:
        return ""


def _items_for(kind: str) -> tuple[str, list[dict], object]:
    """返回 (mode, 条目清单, 路径解析器)。条目含 name/modified/size 与各 kind 的附加字段。"""
    if kind in ("drug_sources", "disease_sources"):
        spec = drug_spec() if kind == "drug_sources" else disease_spec()
        ns = "drug" if kind == "drug_sources" else "disease"
        source_df = list_cache_files(ns)
        merged_df = list_cache_files("merged")
        merged_df = merged_df[merged_df["File name"].str.startswith(f"merged_{ns}_")]
        names = [*source_df["File name"].tolist(), *merged_df["File name"].tolist()]

        def resolver(name: str) -> str:
            return spec_file_path(name, spec)

        items = []
        for n in names:
            path = resolver(n)
            status = merged_status(n, spec)
            items.append({
                "name": n, "modified": _list_mtime(path),
                "size": os.path.getsize(path) if os.path.exists(path) else 0,
                "is_merged": status is not None,
                "merged_status": status,
                "source": next((p.rstrip("_") for p in spec["source_prefixes"]
                                if n.startswith(p)), "merged"),
            })
        return "files", items, resolver

    if kind == "ppi":
        df = list_cache_files("ppi")

        def resolver(name: str) -> str:
            return os.path.join(config.PPI_CACHE_DIR, name)

        items = [{"name": n, "modified": m, "size": os.path.getsize(resolver(n))
                  if os.path.exists(resolver(n)) else 0}
                 for n, m in zip(df["File name"].tolist(), df["Modified time"].tolist())]
        return "files", items, resolver

    if kind == "enrichment":
        df = list_cache_files("enrichment")

        def resolver(name: str) -> str:
            return os.path.join(config.ENRICHMENT_CACHE_DIR, name)

        items = []
        for n, m in zip(df["File name"].tolist(), df["Modified time"].tolist()):
            info = parse_enrichment_filename(n) or {}
            items.append({
                "name": n, "modified": m,
                "size": os.path.getsize(resolver(n)) if os.path.exists(resolver(n)) else 0,
                "display": enrichment_display_name(n),
                "source": info.get("source", ""),
                "side": info.get("side") or "",
                "gsea": info.get("kind") == "gsea",
            })
        return "files", items, resolver

    if kind == "stars":
        star_dir = config.MD_STAR_DIR
        names = sorted(n for n in os.listdir(star_dir)
                       if os.path.isdir(os.path.join(star_dir, n))) \
            if os.path.isdir(star_dir) else []
        return "dirs", [{"name": n, "modified": _list_mtime(os.path.join(star_dir, n))}
                        for n in names], lambda n: os.path.join(star_dir, n)

    if kind == "downloads":
        df = list_download_dirs()
        excluded = {os.path.basename(config.MD_RESULT_DIR),
                    os.path.basename(config.MD_STAR_DIR)}
        df = df[~df["Directory"].isin(excluded)]
        # name 保持 downloads 相对路径(zip/删除的校验键,分析/协同页也在传);
        # display 剥掉 ppi-results/ 聚合层前缀,清单只呈现 run 目录名。
        items = [{"name": d, "modified": m,
                  "display": d[len("ppi-results/"):] if d.startswith("ppi-results/") else d}
                 for d, m in zip(df["Directory"].tolist(),
                                 df["Modified time"].tolist())]
        return "dirs", items, lambda n: os.path.join(config.DOWNLOADS_DIR, n)

    raise HTTPException(status_code=404, detail=f"unknown collection '{kind}'")


@router.get("/api/collections/{kind}")
def collection_items(kind: str) -> dict:
    mode, items, _ = _items_for(kind)
    return {"kind": kind, "mode": mode, "items": items}


class NamesBody(BaseModel):
    names: list[str]


def _validated_paths(kind: str, names: list[str]) -> list[tuple[str, str]]:
    """把选中名解析为 (清单名, 路径),确认都在集合清单内(防越权删除/打包)。
    清单全名与裸名(去 "ppi-results/" 聚合层前缀)都收:管理页传全名,分析页
    拿的是 directory_rel 剥出来的裸名;路径一律按清单条目解析,裸名跨根歧义
    时拒绝。known 只来自 os.listdir,不含分隔符穿越序列,足够防越权。"""
    mode, items, resolver = _items_for(kind)
    full = {it["name"] for it in items}
    by_bare: dict[str, str | None] = {}
    for it in items:
        b = os.path.basename(it["name"])
        if b not in by_bare:
            by_bare[b] = it["name"]
        elif by_bare[b] != it["name"]:
            by_bare[b] = None          # 同名条目出现在两处,裸名不可消歧
    out = []
    for n in names:
        name = n if n in full else by_bare.get(n)
        if not name:
            raise HTTPException(status_code=400, detail=f"invalid name: {n!r}")
        out.append((name, resolver(name)))
    return out


@router.post("/api/collections/{kind}/zip")
def collection_zip(kind: str, body: NamesBody) -> Response:
    mode, _, resolver = _items_for(kind)
    pairs = _validated_paths(kind, body.names)
    if not pairs:
        raise HTTPException(status_code=400, detail="empty selection")
    if mode == "files":
        data = make_ppi_zip({n: p for n, p in pairs})
    else:
        data = make_dirs_zip({n: p for n, p in pairs})
    return _zip_response(data, _ZIP_NAMES[kind])


@router.post("/api/collections/{kind}/delete")
def collection_delete(kind: str, body: NamesBody) -> dict:
    mode, _, resolver = _items_for(kind)
    pairs = _validated_paths(kind, body.names)
    deleted, failed = 0, []
    for n, path in pairs:
        try:
            if mode == "dirs":
                shutil.rmtree(path, ignore_errors=True)
                if os.path.isdir(path):  # rmtree 静默失败(Windows 目录被占用)按失败处理
                    raise OSError(path)
            else:
                os.remove(path)
        except OSError:
            failed.append(n)
        else:
            deleted += 1
    return {"deleted": deleted, "failed": failed}
