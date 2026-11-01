/** NPAP API 的类型化访问层:与 WY/api.py 的端点一一对应。 */

export interface LogFile {
  name: string
  modified: string
  size: number
  is_current: boolean
}

export interface LogContent {
  name: string
  content: string
  truncated: boolean
  size: number
}

export type CollectionKind =
  | 'drug_sources' | 'disease_sources'
  | 'ppi' | 'enrichment' | 'stars' | 'downloads'

export interface CollectionItem {
  name: string
  modified: string
  size: number
  /* 管理页(drug/disease_sources)附加 */
  is_merged?: boolean
  source?: string
  /* 富集页附加 */
  display?: string
  side?: string
  gsea?: boolean
}

export interface CollectionListing {
  kind: CollectionKind
  mode: 'files' | 'dirs'
  items: CollectionItem[]
}

export interface DeleteResult {
  deleted: number
  failed: string[]
}

export interface JobState {
  job_id: string
  status: 'running' | 'done' | 'error'
  description: string
  progress: { fraction: number | null; text: string; [k: string]: unknown }
  result: Record<string, unknown> | null
  error: string | null
  trace?: string
}

export type Row = Record<string, unknown>

/** 网络邻近度检验的落网覆盖(LCC 内命中数 / 集合总数)。 */
export interface ProximityCoverage {
  drug_in_lcc?: number
  drug_total?: number
  disease_in_lcc?: number
  disease_total?: number
}

/** 网络邻近度检验(core.network_algorithms;not_built / 异常哨兵只含部分字段,
 * 故除 ok/error 外均为可选;z / p_excl 在零分布退化时为 null)。 */
export interface ProximityResult {
  ok: boolean
  error: string | null
  d_c?: number | null
  z?: number | null
  p_empirical?: number | null
  d_c_excl?: number | null
  z_excl?: number | null
  p_excl?: number | null
  n_excl?: number
  n_perm?: number
  seed?: number
  mean_random?: number | null
  std_random?: number | null
  mean_random_excl?: number
  lcc_nodes?: number
  cached?: boolean
  coverage?: ProximityCoverage
}

/** RWR 单行(基因得分;annotate 启用时非种子行附带 annotated 标注)。 */
export interface RwrRow {
  gene: string
  score: number
  annotated?: boolean
}

/** RWR 度匹配校正的富集检验;零分布退化(degenerate)时 z / p 为 null。 */
export interface RwrEnrichment {
  n_annotated: number
  obs_mean: number
  null_mean?: number | null
  n_perms: number
  degenerate: boolean
  z: number | null
  p: number | null
}

/** 单侧 RWR 结果(疾病种子 / 药物种子;not_built 哨兵只含部分字段)。 */
export interface RwrSideResult {
  ok: boolean
  error: string | null
  alpha?: number
  n_iter?: number | null
  converged?: boolean | null
  n_seeds?: number | null
  lcc_nodes?: number
  cached?: boolean
  coverage?: { seed_in_lcc?: number; seed_total?: number }
  enrichment?: RwrEnrichment | null
  seeds?: RwrRow[]
  top_nonseed?: RwrRow[]
}

export interface PpiResult {
  prefix: string
  directory_rel: string
  params: (string | number)[]
  hub_rule: string
  hub_method_eff: string
  genes_count: number
  nodes: number
  edges: number
  hubs: number
  expansion: { seeds: number; expanded: number } | null
  hub_table: Row[]
  views: string[]
  artifacts: Record<string, string>
  enrichment: { source: string; reason: string; rows: Row[] }
  gsea: Record<string, { reason: string; stats: Row; rows: Row[] }>
  proximity: ProximityResult | null
  rwr: Record<string, RwrSideResult> | null
  modules: {
    summary: Row[]
    failures: unknown[]
    no_sig: number[]
    /** 单基因模块数(汇总表只含 ≥2 基因的模块,页面提示"已隐藏 N 个") */
    hidden_single?: number
    /** 模块富集长表(按 module_id, p_value 排序,截断前 100 行)与总量统计 */
    enr_rows?: Row[]
    enr_total?: number
    enr_modules?: number
  } | null
}

export interface PpiAnalyzeBody {
  mode: 'names' | 'genes'
  drug: string
  disease: string
  genes_text: string
  score: number
  network_type: string
  rwr_expand: number
  hub_top_n: number
  hub_method: string
}

export interface InteractomeStatus {
  built: boolean
  version?: string
  nodes?: number
  edges?: number
  confidence?: number
  built_at?: string
  age?: number | null
  stale?: boolean
}

/* ---- GEO 表达谱 + ML(阶段 3c)共享类型:字段与 core 一一对应 ---- */

/** 数据集可用性判定(core.expression.assess_expression_dataset)。 */
export interface MlAssess {
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

/** /api/ml/overview 的 job result。 */
export interface MlOverviewResult {
  overview: Row[]
  counts_rows?: Row[]
  assess?: MlAssess
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
export interface MlStaleReason {
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

/* ---- 数据迁移(设置页:结果导出 / 两段式导入) ---- */

export type TransferScope = 'analysis' | 'star' | 'geo' | 'ui_settings'

/** 导入预览的单条目(kind 为目录类范围时含文件数/体积/冲突标记)。 */
export interface TransferItemReport {
  kind: Exclude<TransferScope, 'ui_settings'>
  name: string
  files: number
  bytes: number
  conflict: boolean
  invalid?: boolean
  invalid_reason?: string
}

/** 上传解包后的预览报告;warnings 为导出端与本机算法/数据版本差异提示。 */
export interface TransferImportReport {
  manifest: { app_version: string; created_at: string; schema: Record<string, unknown> }
  items: TransferItemReport[]
  ui_settings: { present: boolean; keys?: number; conflict?: boolean }
  warnings: string[]
}

/** commit 的落盘报告;失败项留在会话里可调整决策重试。 */
export interface TransferCommitReport {
  applied: { kind: string; name: string; keys?: string[] }[]
  skipped: { kind: string; name: string }[]
  failed: { kind: string; name: string; reason: string }[]
}

export interface TransferDecisionBody {
  kind: string
  name: string
  action: 'overwrite' | 'skip' | 'rename'
  new_name?: string
}

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) throw new Error(`HTTP ${r.status}`)
  return (await r.json()) as T
}

/** 同 json(),但把后端 400/413 的 detail 文本带出来(导入包校验错误要给用户看)。 */
async function jsonWithDetail<T>(r: Response): Promise<T> {
  if (!r.ok) {
    let detail = ''
    try {
      detail = ((await r.json()) as { detail?: string })?.detail ?? ''
    } catch {
      /* 响应体不是 JSON 时退回状态码 */
    }
    throw new Error(detail || `HTTP ${r.status}`)
  }
  return (await r.json()) as T
}

async function post<T>(url: string, body: unknown): Promise<T> {
  return fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body)
  }).then((r) => json<T>(r))
}

/** 从 Content-Disposition 取服务端文件名:优先 RFC 5987 的 filename*(中文药名场景),
 * 退回普通 filename="…";都缺失返回 null,由调用方兜底命名。 */
function attachmentName(disposition: string | null): string | null {
  if (!disposition) return null
  const star = /filename\*\s*=\s*(?:UTF-8|utf-8)''([^;]+)/.exec(disposition)
  if (star) {
    try {
      return decodeURIComponent(star[1].trim())
    } catch {
      /* 非法百分号编码时退回普通 filename */
    }
  }
  const plain = /filename\s*=\s*"([^"]*)"/.exec(disposition) ??
                /filename\s*=\s*([^;]+)/.exec(disposition)
  return plain ? plain[1].trim() : null
}

export const api = {
  /* ---- Logs(阶段 1) ---- */
  logsFiles: (): Promise<{ files: LogFile[] }> =>
    fetch('/api/logs/files').then((r) => json(r)),

  logsContent: (name: string, maxBytes = 2_000_000): Promise<LogContent> =>
    fetch(`/api/logs/content?name=${encodeURIComponent(name)}&max_bytes=${maxBytes}`).then(
      (r) => json(r)
    ),

  logsLevel: (): Promise<{ level: string }> =>
    fetch('/api/logs/level').then((r) => json(r)),

  setLogsLevel: (level: 'INFO' | 'DEBUG'): Promise<{ level: string }> =>
    post('/api/logs/level', { level }),

  /* ---- 文件集合(阶段 2) ---- */
  collection: (kind: CollectionKind): Promise<CollectionListing> =>
    fetch(`/api/collections/${kind}`).then((r) => json(r)),

  zip: async (kind: CollectionKind, names: string[]): Promise<Blob> => {
    const r = await fetch(`/api/collections/${kind}/zip`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ names })
    })
    if (!r.ok) throw new Error(`HTTP ${r.status}`)
    return r.blob()
  },

  remove: (kind: CollectionKind, names: string[]): Promise<DeleteResult> =>
    post(`/api/collections/${kind}/delete`, { names }),

  /* ---- 来源刷新(管理页) ---- */
  refreshSources: (
    kind: 'drug' | 'disease',
    filenames: string[] | null,
    hours?: number
  ): Promise<{ job_id: string }> =>
    post(`/api/sources/${kind}/refresh`, { filenames, hours }),

  job: (jobId: string): Promise<JobState> =>
    fetch(`/api/jobs/${jobId}`).then((r) => json(r)),

  /* ---- 背景互作组 ---- */
  interactome: (): Promise<InteractomeStatus> =>
    fetch('/api/interactome').then((r) => json(r)),

  buildInteractome: (): Promise<{ job_id: string }> =>
    post('/api/interactome/build', {}),

  /* ---- PPI 分析(阶段 3a) ---- */
  analyzePpi: (body: PpiAnalyzeBody): Promise<{ job_id: string }> =>
    post('/api/ppi/analyze', body),

  rerunModules: (directoryRel: string, prefix: string): Promise<{ job_id: string }> =>
    post('/api/ppi/rerun-modules', { directory_rel: directoryRel, prefix }),

  rerunGsea: (
    side: 'drug' | 'disease',
    drug: string,
    disease: string,
    prefix: string
  ): Promise<{ job_id: string }> =>
    post('/api/ppi/rerun-gsea', { side, drug, disease, prefix }),

  /* ---- 主页四库检索 + 设置(阶段 3 补全) ---- */
  retrieveTargets: (drug: string, disease: string): Promise<{ job_id: string }> =>
    post('/api/targets/retrieve', { drug, disease }),

  /* 主页检索结果下载:多项后端打 zip,单项直接回对应文件(items 键见 api.py _EXPORT_ITEMS)。
   * 文件名以响应头 Content-Disposition 为准(命名真相源在后端),取不到时用兜底名。 */
  targetsExport: async (jobId: string, items: string[]): Promise<{ blob: Blob; filename: string }> => {
    const r = await fetch(`/api/targets/export?job_id=${encodeURIComponent(jobId)}` +
      `&items=${encodeURIComponent(items.join(','))}`)
    if (!r.ok) throw new Error(`HTTP ${r.status}`)
    const blob = await r.blob()
    return { blob, filename: attachmentName(r.headers.get('content-disposition')) ?? 'targets_download' }
  },

  /* ---- 主页页脚统计(本机缓存实数) ---- */
  homeStats: (): Promise<{
    databases: number
    compounds: number
    targets: number
    pathways: number
  }> => fetch('/api/stats').then((r) => json(r)),

  settings: (): Promise<SettingsSnapshot> => fetch('/api/settings').then((r) => json(r)),

  saveApiEnv: (changes: Record<string, string>): Promise<{ applied: string[] }> =>
    post('/api/settings/api-env', { changes }),

  setStringVersion: (version: string): Promise<{ version: string; changed: boolean }> =>
    post('/api/settings/string-version', { version }),

  saveUiSettings: (changes: Record<string, number | boolean | string | string[]>): Promise<{
    snapshot: SettingsSnapshot
    changed: string[]
    purged_merged: number
  }> => post('/api/settings/ui', { changes }),

  /* ---- 数据迁移(设置页) ---- */
  transferExport: (scopes: TransferScope[]): Promise<{ job_id: string }> =>
    post('/api/transfer/export', { scopes }),

  transferDownload: async (filename: string): Promise<Blob> => {
    const r = await fetch(`/api/transfer/export/${encodeURIComponent(filename)}`)
    if (!r.ok) throw new Error(`HTTP ${r.status}`)
    return r.blob()
  },

  transferImport: (file: File): Promise<{ sid: string; report: TransferImportReport }> =>
    fetch('/api/transfer/import', {
      method: 'POST',
      headers: { 'Content-Type': 'application/zip' },
      body: file
    }).then((r) => jsonWithDetail<{ sid: string; report: TransferImportReport }>(r)),

  transferCommit: (
    sid: string,
    decisions: TransferDecisionBody[]
  ): Promise<{ report: TransferCommitReport; ui_snapshot: SettingsSnapshot | null }> =>
    post(`/api/transfer/import/${encodeURIComponent(sid)}/commit`, { decisions }),

  transferCancel: (sid: string): Promise<{ cancelled: boolean }> =>
    fetch(`/api/transfer/import/${encodeURIComponent(sid)}`, { method: 'DELETE' }).then(
      (r) => json(r)
    ),

  /* ---- 分子对接(阶段 3b) ---- */
  dockingSession: (): Promise<{ sid: string }> =>
    post('/api/docking/session', {}),

  dockingSources: (): Promise<{
    ppi_dirs: { name: string; display: string }[]
    enr_files: { name: string; display: string }[]
  }> => fetch('/api/docking/sources').then((r) => json(r)),

  dockingSource: (
    sid: string,
    body: { mode: 'ppi' | 'enr' | 'manual'; name?: string; manual_genes?: string }
  ): Promise<{ label: string; drug: string; disease: string; hub_count: number;
               enr_count: number; candidates: string[] }> =>
    post(`/api/docking/${sid}/source`,
         { mode: body.mode, name: body.name ?? '', manual_genes: body.manual_genes ?? '' }),

  dockingCheckStructures: (
    sid: string,
    body: { use_struct: boolean; allow_af: boolean }
  ): Promise<{ job_id: string | null; direct?: Record<string, unknown> }> =>
    post(`/api/docking/${sid}/check-structures`, body),

  dockingLigand: (
    sid: string,
    body: { source: 'pubchem' | 'smiles'; query?: string; smiles?: string }
  ): Promise<{ ok: boolean; error?: string; cid?: number | null; lig_pdbqt?: string }> =>
    post(`/api/docking/${sid}/ligand`, body),

  dockingReceptorCandidates: (
    sid: string,
    body: { gene: string; allow_af: boolean }
  ): Promise<{ error?: string; rows: { pdb_id: string; resolution: string; method: string;
             title: string; source: string; flag: string; note: string }[]; af_enabled?: boolean }> =>
    post(`/api/docking/${sid}/receptor/candidates`, body),

  dockingReceptorPrepare: (
    sid: string,
    body: { gene: string; pdb_id: string }
  ): Promise<{ job_id: string }> =>
    post(`/api/docking/${sid}/receptor/prepare`, body),

  dockingBox: (
    sid: string,
    body: { infer: boolean; center?: number[]; size?: number[] }
  ): Promise<{ box: { center: number[]; size: number[] }; clamped: boolean;
               inferred_ligands: string[] | null; residues: string[] }> =>
    post(`/api/docking/${sid}/box`, body),

  dockingRun: (sid: string): Promise<{ job_id: string }> =>
    post(`/api/docking/${sid}/run`, {}),

  dockingSave: (
    sid: string,
    body: { pose_mode: number; conclusion: string; drug: string; disease: string }
  ): Promise<{ result_dir: string }> =>
    post(`/api/docking/${sid}/save`, body),

  dockingStar: (sid: string, on: boolean): Promise<{ starred: boolean }> =>
    post(`/api/docking/${sid}/star`, { on }),

  dockingInteractions: (sid: string, pose: number): Promise<{ rows: Row[] }> =>
    fetch(`/api/docking/${sid}/interactions?pose=${pose}`).then((r) => json(r)),

  dockingView3dUrl: (sid: string, pose: number, cutoff: number, focus: boolean): string =>
    `/api/docking/${sid}/view3d?pose=${pose}&cutoff=${cutoff}&focus=${focus ? 1 : 0}`,

  dockingSave3d: (
    sid: string,
    body: { pose: number; caption: string }
  ): Promise<{ path: string; artifact_rel: string; offline_ok: boolean }> =>
    post(`/api/docking/${sid}/save3d`, body),

  /* ---- GEO 表达谱 + ML 共识筛选(阶段 3c) ---- */
  mlDiseases: (): Promise<{ diseases: string[]; drugs: string[] }> =>
    fetch('/api/ml/diseases').then((r) => json(r)),

  mlSearch: (disease: string, limit: number): Promise<{ results: Row[]; disease: string }> =>
    post('/api/ml/search', { disease, limit }),

  mlProbe: (gse: string, platformInfo: { platform: string; type: string }):
    Promise<MlProbeResult> =>
    post('/api/ml/probe', { gse, platform: platformInfo.platform, type_: platformInfo.type }),

  mlOverview: (gse: string): Promise<{ job_id: string }> =>
    post('/api/ml/overview', { gse }),

  mlRanked: (body: { disease: string; drug: string }): Promise<{
    files: { label: string; base_disease: string; drug: string }[]
  }> => post('/api/ml/ranked', body),

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

/** 产物直出:downloads 相对路径 → 伺服 URL(download=true 触发附件下载)。 */
export const artifactsUrl = (rel: string, download = false): string =>
  `/api/artifacts?rel=${encodeURIComponent(rel)}${download ? '&download=1' : ''}`
