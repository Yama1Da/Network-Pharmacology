"""主页四库检索(App.py 检索职责的 TS 化)。"""

from fastapi import APIRouter
from fastapi import HTTPException
from fastapi.responses import Response
from pydantic import BaseModel
from WY.core import jobs
from WY.core.files import make_bytes_zip
from WY.core.text_utils import safe_name


router = APIRouter()


class RetrieveBody(BaseModel):
    drug: str
    disease: str


_VENN_CACHE: dict[str, dict] = {}   # job_id → {"png": bytes}(任务线程写入,端点读取)
_VENN_ORDER: list[str] = []


def _render_venn_image(drug_set: set, disease_set: set, drug_name: str, disease_name: str,
                       dark: bool = False, fmt: str = "png") -> bytes:
    """venn2 服务端渲染(与 App.py 主页同款;线程安全:用 Figure API 而非 pyplot 全局状态)。
    图内文字保持英文:matplotlib 默认字体不支持中文,中文会显示成方框。
    dark=True 按界面深色主题(墨黑)渲染:深底浅字,供深色模式下 venn 不再是刺眼白块。
    fmt 支持 png / svg(svg 供论文排版二次编辑,matplotlib 文字默认转路径,通用软件可正常显示)。"""
    import io

    from matplotlib.figure import Figure
    from matplotlib_venn import venn2

    bg = "#0e0e0e" if dark else "white"
    fg = "#eaeae8" if dark else "black"
    fig = Figure(figsize=(6, 4), facecolor=bg)
    ax = fig.add_subplot()
    ax.set_facecolor(bg)
    v = venn2(
        subsets=(len(drug_set - disease_set), len(disease_set - drug_set),
                 len(drug_set & disease_set)),
        set_labels=(f"{drug_name} targets", f"{disease_name} targets"),
        ax=ax,
    )
    if v:
        for patch in list(v.set_labels) + [v.get_label_by_id(i) for i in ("10", "01", "11")]:
            if patch is not None:
                patch.set_color(fg)
    if v and v.get_label_by_id("11"):
        v.get_label_by_id("11").set_fontsize(12)  # 交集区域字号
    ax.set_title(f"{drug_name} vs {disease_name} — Venn Diagram", fontsize=14, color=fg)
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=150, bbox_inches="tight", facecolor=bg)
    return buf.getvalue()


def _target_rows(df) -> list[dict]:
    """merged 靶点表 → 主页列表行(JSON 安全:NaN→None,分数保留 4 位;去掉 source_count 冗余列)。"""
    if df is None or df.empty:
        return []
    cols = [c for c in df.columns if c != "source_count"]
    rows: list[dict] = []
    for rec in df[cols].to_dict("records"):
        out: dict = {}
        for k, v in rec.items():
            if v is None or v != v:  # v != v 即 NaN
                out[k] = None
            elif isinstance(v, float):
                out[k] = round(v, 4)
            else:
                out[k] = v
        rows.append(out)
    return rows


@router.post("/api/targets/retrieve")
def targets_retrieve(body: RetrieveBody) -> dict:
    """四库并发检索(ChEMBL/DGIdb/OpenTargets/DisGeNET),写入来源与 merged 缓存;
    分析页的 targets 阶段将直接命中本次缓存。完成后附带主页 venn 图(内存缓存,按 job_id 伺服)。"""
    if not body.drug.strip() or not body.disease.strip():
        raise HTTPException(status_code=400, detail="drug and disease are required")

    venn_holder: dict = {}

    def task(progress) -> dict:
        from WY.core.targets import get_targets_concurrent

        def cb(done: list, total: int) -> None:
            progress(fraction=(len(done) / total) if total else None,
                     text=f"四库检索中 · 已完成 {len(done)}/{total}", done=done)

        drug_df, disease_df, failed = get_targets_concurrent(
            body.drug.strip(), body.disease.strip(), progress_cb=cb)
        drug_set = (set(drug_df["approvedSymbol"].dropna().astype(str))
                    if drug_df is not None and not drug_df.empty else set())
        disease_set = (set(disease_df["approvedSymbol"].dropna().astype(str))
                       if disease_df is not None and not disease_df.empty else set())
        if drug_set or disease_set:
            venn_holder["png"] = _render_venn_image(drug_set, disease_set,
                                                    body.drug.strip(), body.disease.strip())
            # 深色版按需渲染:留存集合与命名,首次请求 theme=dark 时渲染并缓存
            venn_holder["args"] = (drug_set, disease_set,
                                   body.drug.strip(), body.disease.strip())
        common = sorted(drug_set & disease_set)
        drug_rows = _target_rows(drug_df)[:1000]
        disease_rows = _target_rows(disease_df)[:1000]
        common_df = None
        common_rows: list[dict] = []
        if drug_rows and disease_rows:
            common_df = drug_df.drop(columns=["source_count"]).merge(
                disease_df.drop(columns=["source_count"]),
                on="approvedSymbol", how="inner")
            common_rows = _target_rows(common_df)[:1000]
        # 完整表 CSV 留存,供 /api/targets/export 下载(展示行截断到 1000,导出不截断)
        venn_holder["drug_csv"] = _targets_csv_bytes(drug_df)
        venn_holder["disease_csv"] = _targets_csv_bytes(disease_df)
        venn_holder["common_csv"] = _targets_csv_bytes(common_df)
        return {"drug_count": len(drug_set), "disease_count": len(disease_set),
                "common_count": len(common),
                "common_genes": common[:200],
                "drug_rows": drug_rows, "disease_rows": disease_rows,
                "common_rows": common_rows,
                "failed": list(failed)}

    job_id = jobs.submit(task, f"targets retrieve {body.drug}/{body.disease}")
    _VENN_CACHE[job_id] = venn_holder
    _VENN_ORDER.append(job_id)
    if len(_VENN_ORDER) > 10:  # 只保留最近 10 次检索的 venn,防内存增长
        old = _VENN_ORDER.pop(0)
        _VENN_CACHE.pop(old, None)
    return {"job_id": job_id}


@router.get("/api/targets/venn")
def targets_venn(job_id: str, theme: str = "light") -> Response:
    holder = _VENN_CACHE.get(job_id)
    if not holder or not holder.get("png"):
        raise HTTPException(status_code=404, detail="no venn for this job")
    if theme == "dark":
        if "png_dark" not in holder:
            holder["png_dark"] = _render_venn_image(*holder["args"], dark=True)
        return Response(content=holder["png_dark"], media_type="image/png")
    return Response(content=holder["png"], media_type="image/png")


_EXPORT_ITEMS = {
    "drug_csv": "{d}_drug_targets.csv",
    "disease_csv": "{p}_disease_targets.csv",
    "common_csv": "{d}_{p}_common_targets.csv",
    "venn_png": "venn_{d}_{p}.png",
    "venn_svg": "venn_{d}_{p}.svg",
}


def _attachment(fname: str) -> str:
    """Content-Disposition: attachment 头。文件名含非 ASCII(中文药名等)时补 RFC 5987 的
    filename*:FastAPI 响应头按 latin-1 编码,直接放中文会 500。"""
    if fname.isascii():
        return f'attachment; filename="{fname}"'
    from urllib.parse import quote

    fallback = fname.encode("ascii", "replace").decode().replace('"', "_")
    return f'attachment; filename="{fallback}"; filename*=UTF-8\'\'{quote(fname)}'


def _targets_csv_bytes(df) -> bytes:
    """merged 靶点表 → CSV 字节(UTF-8 BOM,Excel 直开;去 source_count 冗余列)。空表回空字节。"""
    if df is None or df.empty:
        return b""
    return df.drop(columns=["source_count"], errors="ignore").to_csv(
        index=False, float_format="%.4g").encode("utf-8-sig")


@router.get("/api/targets/export")
def targets_export(job_id: str, items: str) -> Response:
    """主页检索结果下载:items 逗号分隔(drug_csv/disease_csv/common_csv/venn_png/venn_svg)。
    单项直接回文件,多项打 zip(core.files.make_bytes_zip 内存打包)。
    CSV 为完整 merged 表(不受主页 1000 行展示截断);venn 取浅色版(论文用途),svg 首次请求时渲染进缓存。"""
    holder = _VENN_CACHE.get(job_id)
    if not holder:
        raise HTTPException(status_code=404, detail="no retrieve result for this job")
    want = [s.strip() for s in items.split(",") if s.strip()]
    if not want or any(w not in _EXPORT_ITEMS for w in want):
        raise HTTPException(status_code=400, detail="items must be a comma list of: "
                            + ", ".join(sorted(_EXPORT_ITEMS)))
    resolved: dict[str, bytes] = {}
    for w in want:
        if w == "venn_png":
            data = holder.get("png")
        elif w == "venn_svg":
            if "svg" not in holder and "args" in holder:
                holder["svg"] = _render_venn_image(*holder["args"], fmt="svg")
            data = holder.get("svg")
        else:
            data = holder.get(w)
        if not data:
            raise HTTPException(status_code=404, detail=f"{w} is empty for this job")
        resolved[w] = data
    if "args" not in holder:
        raise HTTPException(status_code=404, detail="no targets for this job")
    d, p = safe_name(holder["args"][2]), safe_name(holder["args"][3])
    if len(resolved) == 1:
        w, data = next(iter(resolved.items()))
        fname = _EXPORT_ITEMS[w].format(d=d, p=p)
        media = {"png": "image/png", "svg": "image/svg+xml"}.get(
            fname.rsplit(".", 1)[1], "text/csv")
    else:
        fname = f"{d}_{p}_targets.zip"
        data = make_bytes_zip({_EXPORT_ITEMS[w].format(d=d, p=p): b for w, b in resolved.items()})
        media = "application/zip"
    return Response(content=data, media_type=media,
                    headers={"Content-Disposition": _attachment(fname)})
