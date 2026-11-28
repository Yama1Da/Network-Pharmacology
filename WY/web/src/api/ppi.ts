/** PPI 分析(阶段 3a)+ 疾病模块多方法一致性(模块C 面板):类型与端点,
 * 字段与 api/ppi.py、api/modules.py 一一对应。 */

import { json, post } from './http'
import type { InFlightJob, Row } from './types'

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

export const ppiApi = {
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

  /* 在途任务状态:重进分析页恢复进度(与 synergyJob/mlJob 同形态) */
  ppiJob: (): Promise<InFlightJob & { mode?: string; drug?: string; disease?: string }> =>
    fetch('/api/ppi/job').then((r) => json(r))
}
