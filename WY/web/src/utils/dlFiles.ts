/** 分析页下载清单的产物名 → 人类可读描述。
 * 后端列运行目录全量文件(除 analysis_meta.json 内部元数据)做**保底**,
 * 认识的后缀叠加人类可读描述;没映射的新文件不会漏——自动回退显示原始文件名。
 * 【维护提醒】新增产物类型时记得两步:①在下面 DL_PATTERNS 加一行后缀→描述键
 * 映射;②i18n 词典 web.analysis.dl_desc.* 补描述(zh/en 各一条)。不加也不会
 * 漏列,该文件显示原始文件名(全量枚举 + 友好名叠加)。 */

const DL_PATTERNS: { suffix: string; key: string }[] = [
  { suffix: '_nodes.csv', key: 'nodes' },
  { suffix: '_edges.csv', key: 'edges' },
  { suffix: '_network.html', key: 'html' },
  { suffix: '_graph.graphml', key: 'graphml' },
  { suffix: '_enrichment.csv', key: 'enrichment_csv' },
  { suffix: '_network_string.png', key: 'string_png' },
  { suffix: '_network_string.svg', key: 'string_svg' },
  { suffix: '_network_cytoscape.png', key: 'cyto_png' },
  { suffix: '_network_cytoscape.svg', key: 'cyto_svg' },
  { suffix: '_network.png', key: 'network_png' },
  { suffix: '_network.svg', key: 'network_svg' },
  { suffix: '_core_targets.png', key: 'core_png' },
  { suffix: '_core_targets.svg', key: 'core_svg' },
  { suffix: '_enrichment_network.png', key: 'dtp_png' },
  { suffix: '_enrichment_network.svg', key: 'dtp_svg' },
  { suffix: '_enrichment_bar.png', key: 'enrichment_bar_png' },
  { suffix: '_enrichment_bar.svg', key: 'enrichment_bar_svg' },
  { suffix: '_enrichment_dot.png', key: 'enrichment_dot_png' },
  { suffix: '_enrichment_dot.svg', key: 'enrichment_dot_svg' },
  { suffix: '_gsea_drug.png', key: 'gsea_drug_png' },
  { suffix: '_gsea_drug.svg', key: 'gsea_drug_svg' },
  { suffix: '_gsea_disease.png', key: 'gsea_disease_png' },
  { suffix: '_gsea_disease.svg', key: 'gsea_disease_svg' },
  { suffix: '_gsea_drug.csv', key: 'gsea_drug_csv' },
  { suffix: '_gsea_disease.csv', key: 'gsea_disease_csv' },
  { suffix: '_rwr_disease_seeds.csv', key: 'rwr_disease_csv' },
  { suffix: '_rwr_drug_seeds.csv', key: 'rwr_drug_csv' },
  { suffix: '_proximity.csv', key: 'proximity_csv' },
  { suffix: '_query_genes.csv', key: 'query_genes_csv' },
  { suffix: '_modules.csv', key: 'modules_csv' },
  { suffix: '_module_enrichment.csv', key: 'module_enrichment_csv' },
  { suffix: 'reach_compare.png', key: 'reach_compare_png' },
  { suffix: '_consensus_targets.csv', key: 'consensus_targets_csv' },
  { suffix: '_consensus_methods.csv', key: 'consensus_methods_csv' },
  { suffix: '_consensus_significance.csv', key: 'consensus_significance_csv' },
]

/** 文件名 → 描述;t 为 vue-i18n 的取词函数(未映射键回退原始文件名)。 */
export function dlDesc(name: string, t: (key: string) => string): string {
  // 最长后缀优先(顺序无关):_module_enrichment.csv 也以 _enrichment.csv 结尾、
  // _enrichment_network.png 也以 _network.png 结尾,短后缀先匹配会错标(2026-11-28 勘误)
  let hit: { suffix: string; key: string } | null = null
  for (const p of DL_PATTERNS) {
    if (name.endsWith(p.suffix) && (!hit || p.suffix.length > hit.suffix.length)) hit = p
  }
  if (!hit) return name
  // 注意 vue-i18n 的 te() 不认扁平点号键:以 t() 原样返回键名判定缺失,回退文件名
  const ik = `web.analysis.dl_desc.${hit.key}`
  const s = t(ik)
  return s === ik ? name : s
}
