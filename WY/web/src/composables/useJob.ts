import { onScopeDispose, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type JobState } from '../api'

/** 长任务的提交-轮询:run() 提交后每 0.6s 轮询一次。
 * fraction/text 为总进度;extra 携带后端附加字段(如 PPI 分析的 step/state 阶段状态)。
 * 默认完成/失败文案走 i18n(在 setup 上下文中调用)。
 * 组件卸载(路由切页,App.vue 按 v-if 销毁视图)自动取消轮询:循环每轮先查
 * 取消标志再发请求,置位即退出且不触发 onDone(2026-11 P0-E:原先每 0.6s
 * 永续请求后端直到任务结束,onDone 闭包与结果不可回收);running 期间重复
 * run() 直接拒绝(重入守卫,机制上消除双循环互写)。 */
export function useJob(onDone?: (result: Record<string, unknown>) => void) {
  const { t } = useI18n()
  const jobId = ref('')
  const running = ref(false)
  const fraction = ref<number | null>(null)
  const text = ref('')
  const extra = ref<Record<string, unknown>>({})
  const error = ref('')
  const result = ref<Record<string, unknown> | null>(null)

  let cancelled = false
  onScopeDispose(() => {
    cancelled = true
  })

  async function run(
    submit: Promise<{ job_id: string }>,
    opts: { doneText?: string; failText?: string } = {}
  ): Promise<void> {
    if (running.value) {
      return
    }
    cancelled = false
    running.value = true
    error.value = ''
    result.value = null
    try {
      const { job_id } = await submit
      jobId.value = job_id
      for (;;) {
        if (cancelled) {
          // 已卸载:静默退出(不触发 onDone),结果由下次进页重新拉取/重跑
          return
        }
        const state: JobState = await api.job(job_id)
        fraction.value = state.progress.fraction
        text.value = state.progress.text
        extra.value = { ...state.progress }
        if (state.status === 'done') {
          result.value = state.result
          text.value = opts.doneText ?? t('web.job.done')
          onDone?.(state.result ?? {})
          break
        }
        if (state.status === 'error') {
          error.value = state.error ?? t('web.job.fail')
          text.value = opts.failText ?? t('web.job.fail')
          break
        }
        await new Promise((r) => setTimeout(r, 600))
      }
    } catch (e) {
      error.value = e instanceof Error ? e.message : String(e)
      text.value = opts.failText ?? t('web.job.fail')
    } finally {
      running.value = false
    }
  }

  return { running, fraction, text, extra, error, result, jobId, run }
}
