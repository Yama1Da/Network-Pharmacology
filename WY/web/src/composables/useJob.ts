import { ref } from 'vue'
import { api, type JobState } from '../api'

/** 长任务的提交-轮询:run() 提交后每 0.6s 轮询一次,进度与结果落到本 composable。 */
export function useJob(onDone?: (result: Record<string, unknown>) => void) {
  const running = ref(false)
  const fraction = ref<number | null>(null)
  const text = ref('')
  const error = ref('')
  const result = ref<Record<string, unknown> | null>(null)

  async function run(
    submit: Promise<{ job_id: string }>,
    opts: { doneText?: string; failText?: string } = {}
  ): Promise<void> {
    running.value = true
    error.value = ''
    result.value = null
    try {
      const { job_id } = await submit
      for (;;) {
        const state: JobState = await api.job(job_id)
        fraction.value = state.progress.fraction
        text.value = state.progress.text
        if (state.status === 'done') {
          result.value = state.result
          text.value = opts.doneText ?? '完成'
          onDone?.(state.result ?? {})
          break
        }
        if (state.status === 'error') {
          error.value = state.error ?? '未知错误'
          text.value = opts.failText ?? '失败'
          break
        }
        await new Promise((r) => setTimeout(r, 600))
      }
    } catch (e) {
      error.value = e instanceof Error ? e.message : String(e)
      text.value = opts.failText ?? '失败'
    } finally {
      running.value = false
    }
  }

  return { running, fraction, text, error, result, run }
}
