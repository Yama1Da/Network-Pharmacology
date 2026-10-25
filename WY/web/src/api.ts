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
  modules: { summary: Row[]; failures: unknown[]; no_sig: number[] } | null
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

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) throw new Error(`HTTP ${r.status}`)
  return (await r.json()) as T
}

async function post<T>(url: string, body: unknown): Promise<T> {
  return fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body)
  }).then((r) => json<T>(r))
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

  /* ---- 主页页脚统计(本机缓存实数) ---- */
  homeStats: (): Promise<{
    databases: number
    compounds: number
    targets: number
    pathways: number
  }> => fetch('/api/stats').then((r) => json(r)),

  settings: (): Promise<{
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
  }> => fetch('/api/settings').then((r) => json(r)),

  saveApiEnv: (changes: Record<string, string>): Promise<{ applied: string[] }> =>
    post('/api/settings/api-env', { changes }),

  setStringVersion: (version: string): Promise<{ version: string; changed: boolean }> =>
    post('/api/settings/string-version', { version }),

  saveUiSettings: (changes: Record<string, number | boolean | string[]>): Promise<{
    snapshot: {
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
    }
    changed: string[]
    purged_merged: number
  }> => post('/api/settings/ui', { changes }),

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
