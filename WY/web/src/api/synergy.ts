/** 分离度协同预测(二期,契约 synergy.v1.1):契约类型与端点,
 * 字段与 core.synergy 一一对应。 */

import { json, post } from './http'
import type { InFlightJob } from './types'

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

export const synergyApi = {
  synergyLibrary: (): Promise<{ drugs: string[]; diseases: string[] }> =>
    fetch('/api/synergy/library').then((r) => json(r)),

  synergyJob: (): Promise<InFlightJob> =>
    fetch('/api/synergy/job').then((r) => json(r)),

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
  }
}
