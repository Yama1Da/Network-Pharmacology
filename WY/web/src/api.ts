/** NPAP API 的类型化访问层:与 WY/api.py 的端点一一对应。
 * 大域按目录拆在 api/ 子目录(http 原语 / 跨域类型 / 各域端点),此处组装,
 * 并继续以 ../api 作为消费者的唯一导入面。 */

import { attachmentName, json, jsonWithDetail, post } from './api/http'
import { mlApi } from './api/ml'
import type { InFlightJob, JobState, Row } from './api/types'

/* 类型出口:跨域共享类型与已外移的 ML 域类型,旧导入路径不变 */
export type { InFlightJob, JobState, Row } from './api/types'
export type {
  MlHistoryRun, MlOverviewResult, MlProbeBatchItem,
  MlProbeResult, MlRunResult, MlRunRow
} from './api/ml'

export interface LogFile {
  name: string
  modified: string
  size: number
  is_current: boolean
}

interface LogContent {
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
  merged_status?: 'current' | 'legacy' | null
  source?: string
  /* 富集页附加 */
  display?: string
  side?: string
  gsea?: boolean
}

interface CollectionListing {
  kind: CollectionKind
  mode: 'files' | 'dirs'
  items: CollectionItem[]
}

interface DeleteResult {
  deleted: number
  failed: string[]
}

/** 网络邻近度检验的落网覆盖(LCC 内命中数 / 集合总数)。 */
interface ProximityCoverage {
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
interface RwrEnrichment {
  n_annotated: number
  obs_mean: number
  null_mean?: number | null
  n_perms: number
  degenerate: boolean
  z: number | null
  p: number | null
}

/** 单侧 RWR 结果(疾病种子 / 药物种子;not_built 哨兵只含部分字段)。 */
interface RwrSideResult {
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

/* 逐靶点可达性标注(W7,2026-11-25):tier 四档独立,"暴露未知"≠"不可达"
 * (方案 §答辩FAQ:过滤可关、默认只标注);reachability=null = 未挂载/缓存未就绪 */
export interface ReachGene {
  r_value: number | null
  tier: string
  /** 真因(2026-11-25 按真因拆词):no_exposure=暴露未知 / no_activity=活性未知 */
  unknown_reason?: string | null
}

export interface TargetsReachability {
  drug: string
  reachability: {
    schema: string
    pk_norm_version: number
    tier: string
    r_primary: number | null
    unknown_reason?: string | null
    n_targets: number
    n_computed: number
  } | null
  reason?: string
  genes: Record<string, ReachGene>
}

/* ---- 疾病模块多方法一致性(二期方案 §4.4,模块C 面板):字段与 api/modules.py 一一对应 ---- */

/** consensus 表单行;秩/共识分的 null = 未被该方法捞入(非成员) */
export interface ConsensusRow {
  gene: string
  rwr_rank: number | null
  diamond_rank: number | null
  neighbor_rank: number | null
  n_methods: number
  consensus_score: number | null
  consistent: boolean
  is_seed: boolean
}

/** 度匹配零模型置换检验的单方法结果(n_perms>0 时才有) */
export interface ConsensusSig {
  p: number
  n_perms: number
  observed_cohesion: number
  null_mean: number
  null_sd: number
}

export interface ModulesConsensusResult {
  ok: boolean
  /** 数据侧失败稳定枚举:no_seeds/disease_required/no_disease_targets/
   * interactome_not_built/rwr_failed——前端映射 i18n */
  reason?: string
  schema?: string
  params: { n: number; k: number; n_perms: number | null }
  /** 退化态:种子入 LCC 数 ≥ N,三方法未扩张,一致子集=种子本身(填入无意义) */
  degenerate?: boolean
  coverage: {
    seeds_total: number
    seeds_in_lcc: number
    seeds_out_of_lcc_count: number
    seeds_out_of_lcc: string[]
  }
  module_sizes: { rwr: number; diamond: number; neighbor: number }
  consistent_genes: number
  /** 一致基因全量清单(填入分析输入按钮用;行表有 _ROWS_CAP 截断,清单没有) */
  consistent_list: string[]
  rows_total: number
  rows: ConsensusRow[]
  significance: Record<string, ConsensusSig> | null
  files_written: string[]
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
  /** W7 招牌图摘要:gate=false 时前端显示结论行;gate=true 看 artifacts 的 reach_compare_png */
  reach_compare: { gate: boolean; judged: number; total: number } | null
  /** 运行目录全量文件清单(除 analysis_meta.json 内部元数据)——下载区逐文件列出 */
  files: { name: string; rel: string }[]
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
  building?: boolean
  job_id?: string
  progress?: JobState['progress']
}

/* ---- PK-DB 镜像(浓度可达性数据底座,二期):字段与 core.pk_mirror.mirror_status 一一对应 ---- */
export interface PkMirrorStatus {
  built: boolean
  studies?: number
  tsv_files?: number
  failures?: number
  missing?: number
  built_at?: string
  age?: number | null
  stale?: boolean
  ttl_days?: number
  building?: boolean
  job_id?: string
  progress?: JobState['progress']
  /** 归一化缓存摘要(reachability 层供给;null/缺省 = 未构建) */
  normalization?: {
    pk_norm_version: number
    conflicts_total?: number
    conflict_keys?: number
    /** 争议名台账全量:名字 → {winner_sid, claimants}(PK 卡展开视图) */
    conflicts?: Record<string, { winner_sid: string; claimants: string[] }>
  } | null
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

/* ---- 分离度协同预测(二期,契约 synergy.v1.1) ---- */

interface SynergyRunParams {
  top_n: number | null
  sort: string | null
  proximity_n_perm: number | null
}

export interface SynergyProx {
  ok: boolean | null
  error: string | null
  d_c: number | null
  z: number | null
  p: number | null
}

export interface SynergyContract {
  schema: string
  pair: { drug_a: string; drug_b: string }
  disease: { name: string; n_targets: number; source: string } | null
  inputs: {
    coverage: { a_total: number; a_in_lcc: number; b_total: number; b_in_lcc: number; d_total: number; d_in_lcc: number }
    n_overlap_ab: number
    power_limited: boolean
    interactome: Record<string, unknown>
  }
  separation: {
    s_ab: number
    d_aa: number
    d_bb: number
    d_ab: number
    null: { mean: number; sd: number; z: number | null; p_lower: number; n_perm: number }
  }
  proximity: { z_a: SynergyProx | null; z_b: SynergyProx | null; z_ab: SynergyProx | null }
  geometry_label: 'overlap' | 'separated' | 'tie'
  label_basis: string
  ranks: null
  /** v1.2 E6 挂载块;null = 未挂载/归一化缓存未就绪/旧 v1.1 归档 */
  reachability: SynergyReachBlock | null
  meta: { computed_at: string; engine: string; run: SynergyRunParams }
}

/** reachability.v1 单药块摘要(逐靶点明细外移运行目录 CSV,result.json 不含)。 */
export interface SynergyReachDrug {
  schema: string
  drug: string
  substance_sid: string | null
  /** reachable | partial | unknown | unreachable(unknown=暴露未知,独立态,≠不可达) */
  tier: string
  r_primary: number | null
  unknown_reason: string | null
  exposure: {
    cmax_molar_M: number | null
    cmax_interventions: string
    n_cmax_obs: number
    fub_value: number | null
    fub_imputed: boolean | null
    /** v1.1:非 null = 主 sid 无暴露记录,数据取自该同义兄弟 sid(台账兜底) */
    synonym_fallback_sid: string | null
  }
  n_targets: number
  n_computed: number
  detail_file: string | null
}

/** v1.2 E6 块:sort_key = |z_sep| × min(w_a,w_b);z=None → null + note "z_null"。 */
export interface SynergyReachBlock {
  drug_a: SynergyReachDrug
  drug_b: SynergyReachDrug
  high_confidence: boolean
  sort_key: number | null
  sort_key_note: string | null
}

/** job 结果载荷 = 契约 + 归档目录名(相对 downloads 集合,可打包下载)。 */
export interface SynergyJobResult extends SynergyContract {
  directory_rel: string
}

/** 结果管理页的运行行摘要(后端 core.synergy.list_results)。 */
export interface SynergyRunRow {
  run: string
  mtime: string
  ok: boolean
  reason: string | null
  drug_a: string
  drug_b: string
  disease: string
  geometry_label: string
  s_ab: number | null
  power_limited: boolean | null
  computed_at: string
  /** v1.2 E6 摘要(旧档/未挂载为 null,前端显示"—") */
  r_tier_a: string | null
  r_tier_b: string | null
  r_high_conf: boolean | null
  r_sort_key: number | null
  r_sort_key_note: string | null
  /** unknown 档真因(按真因拆词:no_exposure/no_activity;判档/旧档为 null) */
  r_unknown_a: string | null
  r_unknown_b: string | null
}

interface SynergyAnalyzeBody {
  drug_a: string
  drug_b: string
  disease?: string
  top_n?: number | null
  n_perm?: number | null
}

/* ---- 数据迁移(设置页:结果导出 / 两段式导入) ---- */

type TransferScope = 'ppi' | 'synergy' | 'star' | 'geo' | 'ui_settings'

/** 导入预览的单条目(kind 为目录类范围时含文件数/体积/冲突标记;'analysis' 为旧版包遗留 kind)。 */
export interface TransferItemReport {
  kind: Exclude<TransferScope, 'ui_settings'> | 'analysis'
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
interface TransferCommitReport {
  applied: { kind: string; name: string; keys?: string[] }[]
  skipped: { kind: string; name: string }[]
  failed: { kind: string; name: string; reason: string }[]
}

interface TransferDecisionBody {
  kind: string
  name: string
  action: 'overwrite' | 'skip' | 'rename'
  new_name?: string
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

  /* ---- PK-DB 镜像(二期) ---- */
  pkMirror: (): Promise<PkMirrorStatus> =>
    fetch('/api/pk-mirror').then((r) => json(r)),

  buildPkMirror: (force = false): Promise<{ job_id: string }> =>
    post('/api/pk-mirror/build', { force }),

  /* ---- 分离度协同预测(二期) ---- */
  synergyLibrary: (): Promise<{ drugs: string[]; diseases: string[] }> =>
    fetch('/api/synergy/library').then((r) => json(r)),

  synergyJob: (): Promise<InFlightJob> =>
    fetch('/api/synergy/job').then((r) => json(r)),

  ppiJob: (): Promise<InFlightJob & { mode?: string; drug?: string; disease?: string }> =>
    fetch('/api/ppi/job').then((r) => json(r)),

  synergyAnalyze: (body: SynergyAnalyzeBody): Promise<{ job_id: string }> =>
    post('/api/synergy/analyze', body),

  synergyResults: (): Promise<{ runs: SynergyRunRow[] }> =>
    fetch('/api/synergy/results').then((r) => json(r)),

  synergyResultsDetail: (run: string): Promise<SynergyContract> =>
    fetch(`/api/synergy/results/detail?run=${encodeURIComponent(run)}`).then((r) => json(r)),

  synergyResultsDelete: (
    runs: string[]
  ): Promise<{ deleted: string[]; skipped: { name: string; reason: string }[] }> =>
    post('/api/synergy/results/delete', { runs }),

  synergyResultsZip: async (runs: string[]): Promise<Blob> => {
    const r = await fetch('/api/synergy/results/zip', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ runs })
    })
    if (!r.ok) throw new Error(`HTTP ${r.status}`)
    return r.blob()
  },

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

  /* 疾病模块多方法一致性(模块C 面板):jobs 异步,产物落运行目录进下载列表;
   * n_perms 网页默认 0(置换检验全量档留 CLI/基准层,性能红线见 api/modules.py) */
  modulesConsensus: (body: {
    directory_rel: string
    prefix: string
    mode: string
    drug?: string
    disease?: string
    genes_text?: string
    n?: number | null
    k?: number | null
    n_perms?: number | null
  }): Promise<{ job_id: string }> => post('/api/modules/consensus', body),

  /* ---- 主页四库检索 + 设置(阶段 3 补全) ---- */
  retrieveTargets: (drug: string, disease: string): Promise<{ job_id: string }> =>
    post('/api/targets/retrieve', { drug, disease }),

  /* 逐靶点可达性标注(W7):药名 → {gene: {r_value, tier}} + 药级摘要;
   * reachability=null(缓存未就绪/无靶点)→ 前端不渲染标注列 */
  targetsReachability: (drug: string): Promise<TargetsReachability> =>
    fetch(`/api/targets/reachability?drug=${encodeURIComponent(drug)}`).then((r) => json(r)),

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
    manual_dirs: { name: string; n_genes: number; mtime: string; genes: string[] }[]
  }> => fetch('/api/docking/sources').then((r) => json(r)),

  dockingSource: (
    sid: string,
    body: { mode: 'ppi' | 'enr' | 'manual' | 'manual-file'; name?: string; manual_genes?: string; manual_file?: string }
  ): Promise<{ label: string; drug: string; disease: string; hub_count: number;
               enr_count: number; candidates: string[] }> =>
    post(`/api/docking/${sid}/source`,
         { mode: body.mode, name: body.name ?? '', manual_genes: body.manual_genes ?? '',
           manual_file: body.manual_file ?? '' }),

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

  // 自带受体上传:原始请求体直传(仿 transferImport,不引 multipart),后端
  // detail(扩展名/大小/体检硬伤)经 jsonWithDetail 直达 UI
  dockingReceptorUpload: (sid: string, gene: string, file: File): Promise<{ job_id: string }> =>
    fetch(
      `/api/docking/${sid}/receptor/upload?gene=${encodeURIComponent(gene)}&name=${encodeURIComponent(file.name)}`,
      { method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: file }
    ).then((r) => jsonWithDetail<{ job_id: string }>(r)),

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
    body: { pose: number; caption: string; cutoff: number }
  ): Promise<{ path: string; artifact_rel: string; offline_ok: boolean }> =>
    post(`/api/docking/${sid}/save3d`, body),

  /* ---- GEO 表达谱 + ML 共识筛选(阶段 3c):方法在 api/ml.ts,此处展开并入 ---- */
  ...mlApi
}

/** 产物直出:downloads 相对路径 → 伺服 URL(download=true 触发附件下载)。 */
export const artifactsUrl = (rel: string, download = false): string =>
  `/api/artifacts?rel=${encodeURIComponent(rel)}${download ? '&download=1' : ''}`
