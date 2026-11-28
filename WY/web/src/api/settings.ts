/** 设置(阶段 3 补全):快照读取、API 环境变量、STRING 版本与 UI 偏好保存,
 * 端点与 api/settings.py 一一对应。 */

import { json, post } from './http'

/* ---- 设置快照(GET /api/settings 与保存/导入设置后的回填共用) ---- */
export interface SettingsSnapshot {
  keys: { disgenet: boolean; ncbi: boolean }
  vina_dir: string
  obabel_dir: string
  cli_results_dir: string
  string_version: string
  known_versions: string[]
  ot_top_n: number
  disgenet_pages: number
  chembl_blacklist: string[]
  chembl_blacklist_options: string[]
  chembl_pchembl_threshold: number
  dgidb_min_score: number
  ot_min_score: number
  chembl_corroboration_required: boolean
  enrichment_pvalue: number
  louvain_resolution: number
  ppi_confidence: number
  pocket_cutoff: number
  /* 统计稳健性/研究语义参数(面板可调,随 ui_settings.json 迁移) */
  rwr_alpha: number
  gsea_permutation_num: number
  proximity_n_perm: number
  rwr_null_perms: number
  /* 分离度协同预测(二期):页面选择器的初始默认值 */
  synergy_top_n: number
  synergy_n_perm: number
  gsea_set_min_size: number
  gsea_set_max_size: number
  gsea_fdr_screen: number
  enrichment_sources: string[]
  exhaustiveness: number
  docking_cpu: number
  num_modes: number
  /* 随搜索深度放宽的有效对接超时(秒),对接页参数行展示用 */
  docking_timeout: number
  ui_font_scale: number
  ui_theme: string
  ui_lang: 'zh' | 'en'
}

export const settingsApi = {
  settings: (): Promise<SettingsSnapshot> => fetch('/api/settings').then((r) => json(r)),

  saveApiEnv: (changes: Record<string, string>): Promise<{ applied: string[] }> =>
    post('/api/settings/api-env', { changes }),

  setStringVersion: (version: string): Promise<{ version: string; changed: boolean }> =>
    post('/api/settings/string-version', { version }),

  saveUiSettings: (changes: Record<string, number | boolean | string | string[]>): Promise<{
    snapshot: SettingsSnapshot
    changed: string[]
    purged_merged: number
  }> => post('/api/settings/ui', { changes })
}
