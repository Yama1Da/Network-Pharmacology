import { computed, ref, watch, type ComputedRef, type Ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type MlOverviewResult, type MlRunResult } from '../api'
import { dropBadge } from '../utils/probeBadgeCache'
import { useJob } from './useJob'

/** ML 运行编排的依赖注入(壳层 refs/computed,不回引其它 composable 实例):
 * 数据上下文(gse/disease/tab)作三重快照代际守卫的观察对象;分组正则与总览
 * 归② 区(换 GSE 连正则与总览一并清——样本命名规则变了);本地矩阵输入、
 * 排序 CSV 勾选与三个调节器参数构成 run 参数包;输入源/可用性判定由壳层
 * ① 区 composable 派生传入。 */
export interface MlRunDeps {
  gse: Ref<string>
  disease: Ref<string>
  tab: Ref<'search' | 'direct' | 'local'>
  casePat: Ref<string>
  ctrlPat: Ref<string>
  overview: Ref<MlOverviewResult | null>
  exprPath: Ref<string>
  labelsPath: Ref<string>
  localCase: Ref<string>
  rankedSel: Ref<Set<string>>
  topK: Ref<number>
  prefilterN: Ref<number>
  bootstrapN: Ref<number>
  hasGse: ComputedRef<boolean>
  useLocal: ComputedRef<boolean>
  localPartial: ComputedRef<boolean>
  sourceConflict: ComputedRef<boolean>
  unusable: ComputedRef<boolean>
  flash: (msg: string) => void
  loadHistory: () => Promise<void>
}

/** ML ③ 运行区(MlView):run 提交(代际守卫)、七阶段进度显示、GEO 回退提示
 * 与在途 run 恢复。强制重下(forceJob)与徽章失效跨 MlView ① 区边界,留壳层。
 * 在 setup 上下文中调用。 */
export function useMlRun(deps: MlRunDeps) {
  const { t } = useI18n()
  const resultData = ref<MlRunResult | null>(null)
  const fallback = ref<{ gse: string; date: string } | null>(null)
  const runJob = useJob((result) => {
    // 失效钩子(2026-11-28 定稿,原绑总览 onDone):run 才真正下载数据集落盘,
    // 本地真实状态可能与归档推断有异——完成后作废该 GSE 徽章回灰点重测;
    // 守卫保证 onDone 只在数据上下文未变时触发,gse 即本次 run 的数据集
    const doneGse = String(result?.gse ?? '').trim().toUpperCase()
    if (doneGse) dropBadge(doneGse)
    resultData.value = result as unknown as MlRunResult
    const r = resultData.value
    if (r.geo_fallback) {
      fallback.value = { gse: String(r.gse), date: String(r.geo_cache_date ?? '') }
    } else {
      fallback.value = null
    }
    deps.flash(t('web.ml.run_done_flash', { gse: r.gse, s: r.elapsed }))
    void deps.loadHistory()
  })

  /* 切换数据上下文时清掉上一轮的展示,避免④区残留旧结果:
   * - 换 GSE:结果/总览/回退提示清空,分组正则一并清(样本命名规则变了),
   *   但三个调节器参数保留——那是用户的分析偏好,与数据集无关;
   * - 换疾病名 / 换输入方式(tab):只清结果区。 */
  watch(deps.gse, () => {
    resultData.value = null
    fallback.value = null
    deps.overview.value = null
    deps.casePat.value = ''
    deps.ctrlPat.value = ''
  })
  watch(deps.disease, () => {
    resultData.value = null
    fallback.value = null
  })
  watch(deps.tab, () => {
    resultData.value = null
    fallback.value = null
  })

  const canRun = computed(() =>
    ((deps.hasGse.value && deps.casePat.value.trim() && deps.ctrlPat.value.trim()) || deps.useLocal.value)
    && !deps.sourceConflict.value && !deps.unusable.value)

  async function run(): Promise<void> {
    if (deps.sourceConflict.value) { deps.flash(t('web.ml.conflict_flash')); return }
    if (deps.localPartial.value) { deps.flash(t('web.ml.partial_flash')); return }
    // 代际守卫(P0-F 及其备案扩展):发起时捕获 GSE/疾病名/输入方式三重快照,
    // 任一变化即视为数据上下文已变,在途结果作废不落④区(watch 只清展示,
    // 不阻断 onDone 回填——守卫补上这段)
    const g0 = deps.gse.value.trim().toUpperCase()
    const d0 = deps.disease.value.trim().toLowerCase()
    const tab0 = deps.tab.value
    await runJob.run(api.mlRun({
      mode: deps.useLocal.value ? 'local' : 'gse',
      gse: g0,
      case_pattern: deps.casePat.value, control_pattern: deps.ctrlPat.value,
      expr: deps.exprPath.value, labels: deps.labelsPath.value, case_label: deps.localCase.value,
      ranked_labels: [...deps.rankedSel.value],
      top_k: deps.topK.value, prefilter_top: deps.prefilterN.value, bootstrap_rounds: deps.bootstrapN.value
    }), {
      doneText: t('web.ml.run_done_text'),
      stale: () => deps.gse.value.trim().toUpperCase() !== g0
        || deps.disease.value.trim().toLowerCase() !== d0
        || deps.tab.value !== tab0
    })
  }

  /* 阶段顺序与标签(标签走 i18n) */
  const STAGE_ORDER = ['fetch', 'group', 'preprocess', 'deg', 'ml', 'shortlist', 'save'] as const
  const STATE_ICON: Record<string, string> = { pending: '○', running: '◐', done: '●' }
  const stageLabels = computed<Record<string, string>>(() => ({
    fetch: t('web.ml.stage_fetch'),
    group: t('web.ml.stage_group'),
    preprocess: t('web.ml.stage_preprocess'),
    deg: t('web.ml.stage_deg'),
    ml: t('web.ml.stage_ml'),
    shortlist: t('web.ml.stage_shortlist'),
    save: t('web.ml.stage_save')
  }))
  const stageStates = computed<Record<string, string>>(() => {
    const map: Record<string, string> = {}
    for (const k of STAGE_ORDER) map[k] = 'pending'
    const p = runJob.extra.value as { step?: string; state?: string }
    if (p?.step && p?.state && p.step in map) map[p.step] = p.state
    return map
  })

  /* 在途 run 恢复(跨页离开任务照跑,重进找回进度与结果——同协同预测语义):
   * 挂载时表单是全新默认态,先采纳任务的 GSE 作为数据上下文(gse 模式),
   * 否则结果回填会撞代际守卫;再按捕获快照续轮询。后端在"完成未超回填窗口"
   * 时也吐 job_id(done)——resume 第一拍命中即回填 ④ 区,同当场跑完体验;
   * local 模式的本地文件路径不可采纳,只续轮询——结果照常落 ④ 区。 */
  async function resumeRunJob(): Promise<void> {
    try {
      const st = await api.mlJob()
      if (!st.job_id || (st.status !== 'running' && st.status !== 'done')) return
      if (st.mode === 'gse' && st.gse) deps.gse.value = st.gse
      const g0 = deps.gse.value.trim().toUpperCase()
      const d0 = deps.disease.value.trim().toLowerCase()
      const tab0 = deps.tab.value
      await runJob.resume(st.job_id, {
        doneText: t('web.ml.run_done_text'),
        stale: () => deps.gse.value.trim().toUpperCase() !== g0
          || deps.disease.value.trim().toLowerCase() !== d0
          || deps.tab.value !== tab0
      })
    } catch { /* 状态不可达:照常空表单 */ }
  }

  return {
    resultData, fallback, runJob, canRun, run, resumeRunJob,
    STAGE_ORDER, STATE_ICON, stageLabels, stageStates
  }
}
