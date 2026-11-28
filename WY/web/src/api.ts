/** NPAP API 的类型化访问层:与 WY/api.py 的端点一一对应。
 * 大域按目录拆在 api/ 子目录(http 原语 / 跨域类型 / 各域端点),此处组装,
 * 并继续以 ../api 作为消费者的唯一导入面。 */

import { json, jsonWithDetail, post } from './api/http'
import { collectionsApi } from './api/collections'
import { dockingApi } from './api/docking'
import { mlApi } from './api/ml'
import { ppiApi } from './api/ppi'
import { retrieveApi } from './api/retrieve'
import { settingsApi } from './api/settings'
import type { InFlightJob, JobState } from './api/types'
import type { SettingsSnapshot } from './api/settings'

/* 类型出口:跨域共享类型与已外移域的类型,旧导入路径不变 */
export type { InFlightJob, JobState, Row } from './api/types'
export type { CollectionItem, CollectionKind } from './api/collections'
export type {
  MlHistoryRun, MlOverviewResult, MlProbeBatchItem,
  MlProbeResult, MlRunResult, MlRunRow
} from './api/ml'
export type {
  ConsensusRow, ConsensusSig, ModulesConsensusResult,
  PpiAnalyzeBody, PpiResult, ProximityResult, RwrRow
} from './api/ppi'
export type { ReachGene, TargetsReachability } from './api/retrieve'
export type { SettingsSnapshot } from './api/settings'

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

  /* ---- 文件集合(阶段 2)+ 来源刷新(管理页):方法在 api/collections.ts,此处展开并入 ---- */
  ...collectionsApi,

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

  /* ---- PPI 分析(阶段 3a)+ 模块一致性:方法在 api/ppi.ts,此处展开并入 ---- */
  ...ppiApi,

  /* ---- 主页四库检索(阶段 3 补全):方法在 api/retrieve.ts,此处展开并入 ---- */
  ...retrieveApi,

  /* ---- 主页页脚统计(本机缓存实数) ---- */
  homeStats: (): Promise<{
    databases: number
    compounds: number
    targets: number
    pathways: number
  }> => fetch('/api/stats').then((r) => json(r)),

  /* ---- 设置(阶段 3 补全):方法在 api/settings.ts,此处展开并入 ---- */
  ...settingsApi,

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

  /* ---- 分子对接(阶段 3b):方法在 api/docking.ts,此处展开并入 ---- */
  ...dockingApi,

  /* ---- GEO 表达谱 + ML 共识筛选(阶段 3c):方法在 api/ml.ts,此处展开并入 ---- */
  ...mlApi
}

/** 产物直出:downloads 相对路径 → 伺服 URL(download=true 触发附件下载)。 */
export const artifactsUrl = (rel: string, download = false): string =>
  `/api/artifacts?rel=${encodeURIComponent(rel)}${download ? '&download=1' : ''}`
