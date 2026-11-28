/** GEO 表达谱 + ML 共识筛选(阶段 3c):类型与端点方法,字段与 core 一一对应。 */

import { json, post } from './http'
import type { InFlightJob, Row } from './types'

/** 数据集可用性判定(core.expression.assess_expression_dataset)。 */
interface MlAssess {
  matrix_empty: boolean
  symbol_level: boolean | null
  shape?: [number, number] | null
  gpl: string
  annot: string | null
  group_fields: string[]
  counts_available: boolean | null
  n_samples: number | null
  usable: boolean
}

/** /api/ml/probe 返回;仅前端请求本身失败时只有 err 字段。 */
export interface MlProbeResult {
  gse?: string
  r?: { size_bytes: number; speed_bps: number; reliable?: boolean } | null
  err?: string
  assess?: MlAssess | null
  annot_missing?: boolean
  gpl?: string
  annot?: string
}

/** /api/ml/probe-batch 的单行结果:徽章原料 + 四态归纳(status 由后端归纳,前端直接消费)。 */
export interface MlProbeBatchItem {
  gse: string
  err: string
  assess: MlAssess | null
  counts_confirmed?: boolean | null
  status: 'ok' | 'unusable' | 'unknown' | 'err'
}

/** 分组选择器字段(总览 job 的 group_picker 载荷):取值含整行命中数 n
 * (与 assign_groups 同口径),overlaps 为子串包含型重叠对(应标红禁选)。 */
interface MlGroupValue {
  value: string
  n: number
}

export interface MlGroupField {
  name: string
  label: string
  values: MlGroupValue[]
  overlaps: { a: string; b: string; n: number }[]
  suggest_case?: string | null
  suggest_control?: string | null
}

interface MlGroupPicker {
  fields: MlGroupField[]
}

/** counts 模式分组选择器(发动机不同的原因见 core.expression.counts_picker_fields):
 * 正则匹配对象是 counts 列名 → 建议物是「取值对应列组的最长公共前缀」而非取值本身。 */
interface MlCountsGroup {
  value: string
  columns: string[]
  suggest: string
}

export interface MlCountsField {
  name: string
  label: string
  groups: MlCountsGroup[]
  suggest_case?: string | null
  suggest_control?: string | null
}

interface MlCountsPicker {
  fields: MlCountsField[]
}

/** /api/ml/overview 的 job result。 */
export interface MlOverviewResult {
  overview: Row[]
  counts_rows?: Row[]
  assess?: MlAssess
  group_picker?: MlGroupPicker
  counts_picker?: MlCountsPicker | null
  counts_match?: { matched: number; total: number } | null
  issues?: string[]
  counts_error?: string
}

/** /api/ml/run 的 job result。 */
export interface MlRunResult {
  gse: string
  tag?: string
  out_dir?: string
  out_dir_rel?: string
  artifacts?: Record<string, string>
  roc_png_rel?: string
  consensus_curve_rel?: string
  elapsed?: number
  meta?: {
    n_samples?: number
    n_genes?: number
    cv_auc?: number | null
    cv_mean?: number | null
    cv_sd?: number | null
    cv_repeats?: number | null
    algo_spearman_mean?: number | null
    [k: string]: unknown
  }
  geo_fallback?: boolean
  geo_cache_date?: string
  funnel?: Row[]
  consensus?: Row[]
  consensus_total?: number
  shortlist?: Row[]
  shortlist_count?: number | null
  deg?: Row[]
}

/** run 陈旧原因(core.ml_pipeline.stale_reasons)。 */
interface MlStaleReason {
  type: string
  old?: string | null
  new?: string | null
  name?: string
}

/** /api/ml/history 的单行(近 3 次摘要)。 */
export interface MlHistoryRun {
  name: string
  auc?: number | null
  n_samples?: number | null
  shortlist_count?: number | null
  ml_version?: string | null
  stale?: MlStaleReason[]
}

/** /api/ml/runs 的单行(结果管理清单)。 */
export interface MlRunRow {
  name: string
  path: string
  dir_rel?: string
  size_bytes: number
  mtime: number
  gse?: string | null
  source?: string | null
  n_samples?: number | null
  auc?: number | null
  shortlist_count?: number | null
  ml_version?: string | null
  stale?: MlStaleReason[]
}

export const mlApi = {
  mlDiseases: (): Promise<{ diseases: string[]; drugs: string[] }> =>
    fetch('/api/ml/diseases').then((r) => json(r)),

  mlSearch: (disease: string, limit: number): Promise<{ results: Row[]; disease: string }> =>
    post('/api/ml/search', { disease, limit }),

  mlProbe: (gse: string, platformInfo: { platform: string; type: string }):
    Promise<MlProbeResult> =>
    post('/api/ml/probe', { gse, platform: platformInfo.platform, type_: platformInfo.type }),

  mlProbeBatch: (rows: { gse: string; platform: string; type: string }[]):
    Promise<{ results: MlProbeBatchItem[] }> =>
    post('/api/ml/probe-batch', {
      rows: rows.map((r) => ({ gse: r.gse, platform: r.platform, type_: r.type }))
    }),

  mlOverview: (gse: string): Promise<{ job_id: string }> =>
    post('/api/ml/overview', { gse }),

  mlRanked: (body: { disease: string; drug: string }): Promise<{
    files: { label: string; base_disease: string; drug: string }[]
  }> => post('/api/ml/ranked', body),

  mlJob: (): Promise<InFlightJob & { gse?: string; mode?: string }> =>
    fetch('/api/ml/job').then((r) => json(r)),

  mlRun: (body: Record<string, unknown>): Promise<{ job_id: string }> =>
    post('/api/ml/run', body),

  mlHistory: (): Promise<{ runs: MlHistoryRun[] }> =>
    fetch('/api/ml/history').then((r) => json(r)),

  mlRuns: (): Promise<{ runs: MlRunRow[] }> =>
    fetch('/api/ml/runs').then((r) => json(r)),

  mlRunsDelete: (body: { paths: string[] }): Promise<{ deleted: string[] }> =>
    post('/api/ml/runs/delete', body),

  mlRunsZip: (paths: string[]): Promise<Blob> =>
    fetch('/api/ml/runs/zip', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ paths })
    }).then((r) => {
      if (!r.ok) throw new Error(`HTTP ${r.status}`)
      return r.blob()
    }),

  mlRunsCleanup: (body: { keep: number }): Promise<{ deleted: number }> =>
    post('/api/ml/runs/cleanup', body),

  mlForceRefresh: (gse: string): Promise<{ job_id: string }> =>
    post('/api/ml/force-refresh', { gse })
}
