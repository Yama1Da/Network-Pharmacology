/** 主页四库检索(阶段 3 补全)与逐靶点可达性标注(W7):端点与 api/retrieve.py 的
 * /api/targets/* 一一对应。 */

import { attachmentName, json, post } from './http'

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

export const retrieveApi = {
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
  }
}
