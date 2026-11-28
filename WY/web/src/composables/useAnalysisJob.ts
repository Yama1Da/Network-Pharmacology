import { computed, ref, type ComputedRef, type Ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type PpiAnalyzeBody, type PpiResult, type Row } from '../api'
import { useJob } from './useJob'
import { useReachAnnotation } from './useReach'

/** 分析页主作业的依赖注入:表单 9 字段(壳层持有——一致性面板 fill 联动与
 * datalist 模板都落在壳层)+ 参数快照(paramsStale 比对基准)+ 页顶 flash。 */
export interface AnalysisJobDeps {
  mode: Ref<'names' | 'genes'>
  drug: Ref<string>
  disease: Ref<string>
  genesText: Ref<string>
  score: Ref<number>
  networkType: Ref<string>
  rwrExpand: Ref<number>
  hubTopN: Ref<number>
  hubMethod: Ref<string>
  snapshot: ComputedRef<(string | number)[]>
  flash: (msg: string) => void
}

/** 分析页主作业与结果数据源:ppi 提交/轮询/恢复、七阶段进度显示、hub 表
 * 可达性标注列与过滤(W7)、paramsStale 参数过期标记。ppi 是全页唯一结果
 * 对象——富集/重跑合并/下载/一致性胶水均从壳层解构的 ppi 读取。
 * 在 setup 上下文中调用。 */
export function useAnalysisJob(deps: AnalysisJobDeps) {
  const { t } = useI18n()
  const ppi = ref<PpiResult | null>(null)

  const paramsStale = computed(() => {
    if (!ppi.value) return false
    const a = deps.snapshot.value
    const b = ppi.value.params
    return a.length !== b.length || a.some((v, i) => String(v) !== String(b[i]))
  })

  /* 可达性标注(W7):names 模式才有药物语境,任务完成后懒取,挂在 hub 表;
   * genes 模式粘贴基因清单无暴露语义,不取数不渲染(共享逻辑见 useReach.ts) */
  const {
    reachFilter, load: loadReach, reset: resetReach, hasReach,
    tier: reachTier, word: reachWord, title: reachTitle, matchesFilter
  } = useReachAnnotation()

  const job = useJob((result) => {
    ppi.value = result as unknown as PpiResult
    const p = ppi.value
    if (String(p?.params?.[0]) === 'names') void loadReach(String(p?.params?.[1] ?? ''))
    else resetReach()
  })

  /* hub 表标注列:行内基因列名固定 gene(metrics_df 口径);过滤默认"全部"保一期 */
  const hubCols = computed<string[]>(() => {
    const base = ppi.value?.hub_table.length ? Object.keys(ppi.value.hub_table[0]) : []
    if (hasReach.value && !base.includes('__reach')) base.push('__reach')
    return base
  })
  const hubRowsShown = computed<Row[]>(() => {
    const rows = ppi.value?.hub_table ?? []
    if (!hasReach.value || reachFilter.value === 'all') return rows
    return rows.filter((r) => matchesFilter(reachTier(r.gene), reachFilter.value))
  })

  const stages = computed<string[]>(() => {
    const base = ['targets', 'string', 'graph', 'enrichment']
    const extra = deps.mode.value === 'names' ? ['gsea', 'network'] : []
    return [...base, ...extra, 'modules', 'save']
  })

  const stageLabels = computed<Record<string, string>>(() => ({
    targets: t('web.analysis.stage_targets'),
    string: t('web.analysis.stage_string'),
    graph: t('web.analysis.stage_graph'),
    enrichment: t('web.analysis.stage_enrichment'),
    gsea: t('web.analysis.stage_gsea'),
    network: t('web.analysis.stage_network'),
    modules: t('web.analysis.stage_modules'),
    save: t('web.analysis.stage_save')
  }))
  const STATE_ICON: Record<string, string> = { pending: '○', running: '◐', done: '●' }

  const stepStates = computed<Record<string, string>>(() => {
    const map: Record<string, string> = {}
    for (const k of stages.value) map[k] = 'pending'
    const p = job.extra.value as { step?: string; state?: string; steps?: Record<string, string> }
    // 后端富集/GSEA/邻近度+RWR/模块检测并发执行,整体上报 steps 映射 → 多阶段同时点亮;
    // steps 缺失时(旧任务/其他作业)退回单 step/state 语义
    if (p?.steps) {
      for (const [k, st] of Object.entries(p.steps)) if (k in map) map[k] = st
    } else if (p?.step && p?.state && p.step in map) {
      map[p.step] = p.state
    }
    return map
  })

  function analyze(): void {
    if (deps.mode.value === 'names' && (!deps.drug.value.trim() || !deps.disease.value.trim())) {
      deps.flash(t('web.analysis.need_names'))
      return
    }
    // genes 模式下 drug/disease 是隐藏表单的残留值,与结果无关:传空,后端 params
    // 快照同样置空,两侧一致,paramsStale 才不会误报参数过期
    const namesMode = deps.mode.value === 'names'
    const body: PpiAnalyzeBody = {
      mode: deps.mode.value, drug: namesMode ? deps.drug.value : '', disease: namesMode ? deps.disease.value : '',
      genes_text: deps.genesText.value,
      score: deps.score.value, network_type: deps.networkType.value, rwr_expand: deps.rwrExpand.value,
      hub_top_n: deps.hubTopN.value, hub_method: deps.hubMethod.value
    }
    void job.run(api.analyzePpi(body), { doneText: t('web.analysis.report_done') })
  }

  /* 在途分析恢复(跨页离开任务照跑,重进找回进度与结果——同协同预测语义):
   * names 模式先采纳任务的分析对象,结果回填后 paramsStale 不会误报参数过期;
   * 后端在"完成未超回填窗口"时也吐 job_id(done)——resume 即回填报告,同当场
   * 跑完体验;genes 模式基因列表不可采纳,只续轮询。 */
  async function resumeAnalyzeJob(): Promise<void> {
    try {
      const st = await api.ppiJob()
      if (!st.job_id || (st.status !== 'running' && st.status !== 'done')) return
      if (st.mode === 'names' && st.drug) {
        deps.mode.value = 'names'
        deps.drug.value = st.drug
        deps.disease.value = st.disease ?? ''
      }
      await job.resume(st.job_id, { doneText: t('web.analysis.report_done') })
    } catch { /* 状态不可达:照常空表单 */ }
  }

  return {
    ppi, paramsStale, job,
    reachFilter, hasReach, reachTier, reachWord, reachTitle,
    hubCols, hubRowsShown,
    stages, stageLabels, STATE_ICON, stepStates,
    analyze, resumeAnalyzeJob
  }
}
