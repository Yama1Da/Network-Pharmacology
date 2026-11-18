"""后端英文文案表(默认语言):与 translations_zh 同键集。TS 前端文案(含顶栏
导航 web.nav.*)全部在 tools/export_i18n_json.py 的 web.* 命名空间;
2026-11-18 死键清理删除 Streamlit 时代遗留的 620 个无引用键(原 640 键),
同日前端导航改用 web.nav.* 后移除 app.btn_*/app.popover_* 复用键。"""

TEXT_EN: dict[str, str] = {
    # Venn 图内文字刻意保持英文，硬编码在 App.py 图表处，相关键已移除。
    # 富集文件友好名称中的 GSEA 侧别标签
    "enrich_file.side_drug": "drug side",
    "enrich_file.side_disease": "disease side",
    "enrich_file.module_n": "module {n}",
    # 管理页共用:全部更新进度条
    "mgmt.prog.sources": "Updating source files ({done}/{total})… last done: {latest}",
    "mgmt.prog.merged": "Rebuilding merged files ({done}/{total})…",
    "mgmt.prog.done": "All files processed",
    "dock.err_invalid_smiles": "Invalid SMILES: {s}",
    "dock.err_rdkit_3d": "RDKit failed to generate 3D conformers",
    "dock.err_mode_range": "mode {mode} out of range ({n} total)",
}
