"""后端中文文案表:仅收后端 t() 在用的错误提示/进度文案。
TS 前端文案(含顶栏导航 web.nav.*)全部在 tools/export_i18n_json.py 的
web.* 命名空间;2026-11-18 死键清理删除 Streamlit 时代遗留的 620 个无引用键
(原 640 键),同日前端导航改用 web.nav.* 后移除 app.btn_*/app.popover_* 复用键。"""

TEXT_ZH: dict[str, str] = {
    # 富集文件友好名称中的 GSEA 侧别标签
    "enrich_file.side_drug": "药物侧",
    "enrich_file.side_disease": "疾病侧",
    "enrich_file.module_n": "模块 {n}",
    # 管理页共用:全部更新进度条
    "mgmt.prog.start": "开始更新来源文件（共 {total} 个）…",
    "mgmt.prog.sources": "正在更新来源文件（{done}/{total}）… 刚完成：{latest}",
    "mgmt.prog.merged": "正在重建整合文件（{done}/{total}）…",
    "mgmt.prog.done": "全部文件已处理",
    "dock.err_invalid_smiles": "无效 SMILES：{s}",
    "dock.err_rdkit_3d": "RDKit 3D 构象生成失败",
    "dock.err_mode_range": "mode {mode} 超出范围（共 {n} 个）",
}
