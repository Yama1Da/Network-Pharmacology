import { computed, ref } from 'vue'
import { api, type JobState } from '../api'

/** 任务中心全局单例:跨页跟踪全部长任务(提交-轮询模型的"看门人")。
 *
 * useJob.run/resume 拿到 job_id 后调用 trackJob() 自动登记,路由取提交时
 * 所在页(location.hash);任务名由后端注册表 description 前缀映射(NAME_KEYS),
 * 尾参(药物/疾病/GSE 等)作为附注展示。本文件只产状态,呈现与交互在
 * TaskCenter.vue:慢轮询(1s,仅对在途任务,全部终态即零请求)、终态事件
 * (组件据此弹气泡)、未读数(白条/角标用)。已在页面内亲眼看完的任务由
 * useJob 的轮询就地 readTask(),不产生未读。 */

export interface CenterTask {
  jobId: string
  route: string                  // 提交时所在页:「查看结果」跳回这里
  description: string            // 后端注册表描述(如 "ppi analyze 黄芪/2型糖尿病")
  nameKey: string                // description 前缀映射的 i18n 键
  nameRest: string               // 前缀后的尾参(药物/疾病/GSE 等),作附注
  status: 'running' | 'done' | 'error'
  fraction: number | null
  text: string
  finishedAt: number | null      // ms;完成时刻(回填窗口倒计时起点)
  error: string | null
  fails: number                  // 连续轮询失败数(404/GC 兜底用)
}

/** 与后端 RESULT_REPLAY_WINDOW_S 对齐:完成后可回填的时长(倒计时展示用)。 */
export const REPLAY_WINDOW_S = 600
const POLL_MS = 1000
const MAX_FAILS = 5             // 连续失败(404=已被 GC / 后端重启):按丢失终态处理

export const tasks = ref<CenterTask[]>([])
export const unread = computed(() => tasks.value.filter(t => t.status !== 'running').length)
/** 最近一次「在途 → 终态」事件;seq 递增,保证同任务重复事件也能触发 watch。 */
export const lastEvent = ref<{ jobId: string; ok: boolean; seq: number } | null>(null)

let seq = 0
let timer = 0

/** description 前缀 → i18n 任务名。长前缀在前(如 receptor upload 先于 receptor)。 */
const NAME_KEYS: ReadonlyArray<readonly [string, string]> = [
  ['ppi analyze', 'web.jobs.ppi'],
  ['ppi rerun modules', 'web.jobs.modules_rerun'],
  ['ppi rerun gsea', 'web.jobs.gsea'],
  ['modules consensus', 'web.jobs.consensus'],
  ['ml overview', 'web.jobs.ml_overview'],
  ['ml run', 'web.jobs.ml_run'],
  ['geo force refresh', 'web.jobs.ml_force'],
  ['docking structure check', 'web.jobs.docking_check'],
  ['docking receptor upload', 'web.jobs.docking_receptor'],
  ['docking receptor', 'web.jobs.docking_receptor'],
  ['docking vina run', 'web.jobs.docking_run'],
  ['pk-db mirror', 'web.jobs.pk'],
  ['interactome build', 'web.jobs.interactome'],
  ['targets retrieve', 'web.jobs.retrieve'],
  ['synergy ', 'web.jobs.synergy'],
  ['drug sources refresh', 'web.jobs.source_refresh'],
  ['disease sources refresh', 'web.jobs.source_refresh'],
  ['results export', 'web.jobs.transfer'],
]

function resolveName(description: string): { key: string; rest: string } {
  const lower = description.toLowerCase()
  for (const [prefix, key] of NAME_KEYS) {
    if (lower.startsWith(prefix)) {
      return { key, rest: description.slice(prefix.length).trim() }
    }
  }
  return { key: 'web.jobs.generic', rest: description }
}

function currentRoute(): string {
  return location.hash.replace(/^#\/?/, '') || 'home'
}

/** 登记任务(useJob.run/resume 拿到 job_id 后调用;同 id 幂等)。 */
export function trackJob(jobId: string): void {
  if (tasks.value.some(t => t.jobId === jobId)) return
  tasks.value.push({
    jobId, route: currentRoute(), description: '',
    nameKey: 'web.jobs.generic', nameRest: '',
    status: 'running', fraction: null, text: '',
    finishedAt: null, error: null, fails: 0,
  })
  if (!timer) timer = window.setInterval(tick, POLL_MS)
}

/** 已读:移除条目(查看结果/查看日志跳页,或页面内亲眼看完时由 useJob 调用)。 */
export function readTask(jobId: string): void {
  tasks.value = tasks.value.filter(t => t.jobId !== jobId)
}

async function tick(): Promise<void> {
  const running = tasks.value.filter(t => t.status === 'running')
  if (running.length === 0) return                    // 全部终态:零请求,循环空转
  await Promise.all(running.map(async (t) => {
    let st: JobState | null = null
    try {
      st = await api.job(t.jobId)
    } catch {
      t.fails += 1                                    // 网络抖动:下一拍再试
      if (t.fails >= MAX_FAILS) {
        // 404(已被 GC)/后端重启:按丢失终态处理,避免永久在途
        t.status = 'error'
        t.error = 'lost'
        t.finishedAt = Date.now()
        lastEvent.value = { jobId: t.jobId, ok: false, seq: ++seq }
      }
      return
    }
    t.fails = 0
    t.description = st.description || t.description
    const name = resolveName(t.description)
    t.nameKey = name.key
    t.nameRest = name.rest
    t.fraction = st.progress.fraction
    t.text = st.progress.text
    if (st.status === 'running') return
    t.status = st.status
    t.error = st.error
    // 后端 finished(秒)优先;缺失则按当前时刻起算(回填窗口倒计时起点)
    t.finishedAt = typeof st.finished === 'number' ? st.finished * 1000 : Date.now()
    lastEvent.value = { jobId: t.jobId, ok: t.status === 'done', seq: ++seq }
  }))
}
