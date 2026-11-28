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
import { synergyApi } from './api/synergy'
import type { JobState } from './api/types'
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
export type {
  SynergyContract, SynergyJobResult, SynergyProx,
  SynergyReachBlock, SynergyReachDrug, SynergyRunRow
} from './api/synergy'

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

  /* ---- 分离度协同预测(二期):方法在 api/synergy.ts,此处展开并入 ---- */
  ...synergyApi,

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
