/** 数据迁移(设置页:结果导出 / 两段式导入),端点与 api/transfer.py 一一对应。 */

import { json, jsonWithDetail, post } from './http'
import type { SettingsSnapshot } from './settings'

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

export const transferApi = {
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
    )
}
