<script setup lang="ts">
/** 任务中心悬浮胶囊:边缘吸附(缝隙)/弹出折叠胶囊/展开面板三态,状态源自
 * useJobCenter 单例(跨页跟踪全部长任务),本组件只做呈现与交互。
 * 交互语义(与 mockups/task-center-options.html 示意稿一致):
 * · 唯一休息形态 = 贴右缘缝隙(12px 探头):空闲灰缝不渲染、运行中黑点沿边框
 *   巡游(SVG dashoffset,动效即运行信号)、有未读黑缝白条(1 条 = 1 个,封顶 5);
 * · 悬停探头带 260ms 意图延迟;点击弹出折叠胶囊;再点展开面板(永远完整在界内);
 * · 完成/失败事件弹出锚定气泡,停 5s 自动缩回,悬停冻结倒计时;
 * · 从面板收起同样停 5s 缩回;拖动沿竖边换停靠位置(localStorage 记忆),向外拖即弹出。
 * 颜色全部走 --c-* 主题变量(深色主题自动反转)。 */
import { computed, nextTick, onBeforeUnmount, onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import {
  tasks, unread, lastEvent, readTask, REPLAY_WINDOW_S, type CenterTask
} from '../composables/useJobCenter'

const { t } = useI18n()

const RETRACT_MS = 5000      // 弹出后停留时长(未读事件与面板收起同规)
const PEEK_DELAY_MS = 260    // 悬停意图延迟(扫过不探头,停住才探头)
const BAR_MAX = 5            // 缝隙白条封顶;真实数量由胶囊文字与面板表达
const DOCK_Y_KEY = 'npap-tc-dock-y'

type Ui = 'dock' | 'pop' | 'panel'
const ui = ref<Ui>('pop')
const dockY = ref(Math.round(window.innerHeight * 0.66))   // 默认中间偏下;localStorage 记忆优先
const dragging = ref(false)
const hovered = ref(false)
const retractAt = ref<number | null>(null)
const bubbleJobId = ref<string | null>(null)
/** 欢迎态:挂载后弹出欢迎词 3s,期间完全被动(不响应悬停/点击/拖动);真实任务事件优先 */
const welcome = ref(true)
const now = ref(Date.now())
const flipBelow = ref(false)

const launcherEl = ref<HTMLElement | null>(null)
const pillEl = ref<HTMLElement | null>(null)
const panelEl = ref<HTMLElement | null>(null)
const pillW = ref(0)
const pillH = ref(0)

const runningCount = computed(() => tasks.value.filter(t => t.status === 'running').length)
const runningAvgPct = computed(() => {
  const list = tasks.value.filter(t => t.status === 'running' && t.fraction != null)
  if (list.length === 0) return null
  return Math.round((list.reduce((s, x) => s + (x.fraction ?? 0), 0) / list.length) * 100)
})
/** 缝隙/胶囊的身份:有未读 > 在跑 > 空闲(空闲也渲染——灰缝是常驻入口) */
const status = computed(() => unread.value > 0 ? 'unread' : runningCount.value > 0 ? 'run' : 'idle')

const pillText = computed(() => {
  if (ui.value === 'panel') return t('web.jobs.collapse')
  if (status.value === 'unread') return t('web.jobs.pending_n', { n: unread.value })
  if (status.value === 'run') {
    const base = t('web.jobs.running_n', { n: runningCount.value })
    return runningAvgPct.value == null ? base : `${base} · ${runningAvgPct.value}%`
  }
  return t('web.jobs.idle')
})

const bubbleTask = computed(() =>
  bubbleJobId.value ? tasks.value.find(x => x.jobId === bubbleJobId.value) ?? null : null)

function remainingSec(task: CenterTask): number {
  const finished = task.finishedAt ?? now.value
  return Math.max(0, Math.floor(REPLAY_WINDOW_S - (now.value - finished) / 1000))
}
function fmt(sec: number): string {
  const s = Math.max(0, sec)
  return `${String(Math.floor(s / 60)).padStart(2, '0')}:${String(s % 60).padStart(2, '0')}`
}
function taskName(task: CenterTask): string {
  return t(task.nameKey) + (task.nameRest ? ` · ${task.nameRest}` : '')
}

/* ---------- 边框巡游:SVG 虚线沿轮廓匀速绕行(仅运行中渲染) ---------- */
function marchPath(w: number, h: number): string {
  // 两条轮廓同为顺时针(缝隙:底边右→左、左侧下→上、顶边左→右),巡游方向一致
  return w <= 13
    ? `M ${w - 0.75} ${h - 0.75} L 8 ${h - 0.75} Q 0.75 ${h - 0.75} 0.75 ${h - 8} L 0.75 8 Q 0.75 0.75 8 0.75 L ${w - 0.75} 0.75`
    : `M 0.75 0.75 L ${w - 0.75} 0.75 L ${w - 0.75} ${h - 0.75} L 0.75 ${h - 0.75} Z`
}
const sliverPath = computed(() => marchPath(12, Math.max(pillH.value, 34)))
const ringPath = computed(() => marchPath(Math.max(pillW.value, 120), Math.max(pillH.value, 34)))
let ro: ResizeObserver | null = null

/* ---------- 拖动与点击切换都在 pointer 事件里完成(5px 位移阈值区分):
 * setPointerCapture 会把合成 click 重定向到捕获元素(launcher 根),挂在子元素
 * 上的 @click 收不到——切换必须留在 pointerup。面板/气泡内的按下不入 press
 * (onDown 早退),对应 pointerup 按"没有按下"忽略,点「查看结果」不会折叠面板。 ---------- */
let press: { x: number; y: number; y0: number; moved: boolean } | null = null
function onDown(e: PointerEvent): void {
  // 欢迎态整段被动;面板/气泡内部的按下与拖动/切换无关(点「查看结果」等不做手势)
  if (welcome.value) return
  if (e.target instanceof Element && e.target.closest('.tc-panel, .tc-bubble')) return
  press = { x: e.clientX, y: e.clientY, y0: dockY.value, moved: false }
  ;(e.currentTarget as HTMLElement).setPointerCapture(e.pointerId)
}
function onMove(e: PointerEvent): void {
  if (!press) return
  const dx = e.clientX - press.x
  const dy = e.clientY - press.y
  if (!press.moved && Math.hypot(dx, dy) < 5) return
  press.moved = true
  dragging.value = true
  if (ui.value === 'dock' && dx < -18) enterPop(false)   // 向外拖 = 弹出
  const max = window.innerHeight - (pillH.value || 34) - 6
  dockY.value = Math.max(6, Math.min(max, press.y0 + dy))
}
function onUp(): void {
  // 没有对应按下(面板/气泡内的点击,pointerdown 已早退):不是手势,直接忽略
  if (!press) return
  const wasDrag = press.moved
  press = null
  dragging.value = false
  if (wasDrag) {
    persistDockY()
    return
  }
  if (ui.value === 'dock') {
    enterPop(true)                                       // 点缝隙 → 弹出(悬停期间不缩)
  } else if (ui.value === 'pop') {
    ui.value = 'panel'                                   // 再点 → 展开面板
    bubbleJobId.value = null
  } else {
    ui.value = 'pop'                                     // 面板收起 → 停 5s 缩回
    retractAt.value = Date.now() + RETRACT_MS
  }
}
function enterPop(withCountdown: boolean): void {
  ui.value = 'pop'
  retractAt.value = withCountdown ? null : Date.now() + RETRACT_MS
}
function onEnter(): void {
  if (welcome.value) return                              // 欢迎词倒计时不冻结
  retractAt.value = null
}
function onLeave(): void {
  if (welcome.value) return
  if (ui.value === 'pop') retractAt.value = Date.now() + RETRACT_MS
}

function collapse(): void {
  ui.value = 'pop'
  retractAt.value = Date.now() + RETRACT_MS
}

/** 「查看结果/查看日志」:跳回提交页(对应视图的 resume/装载逻辑负责回填),并已读 */
function openTask(task: CenterTask): void {
  location.hash = `#/${task.route}`
  readTask(task.jobId)
}

function persistDockY(): void {
  try { localStorage.setItem(DOCK_Y_KEY, String(dockY.value)) } catch { /* 私有模式等:放弃记忆 */ }
}
function clampDockY(): void {
  const max = window.innerHeight - (pillH.value || 34) - 6
  dockY.value = Math.max(6, Math.min(max, dockY.value))
}

/* ---------- 心跳:缩回计时 + 倒计时刷新(就地更新,不重排面板 DOM) ---------- */
let timer = 0
let welcomeTimer = 0
onMounted(() => {
  try {
    const raw = parseFloat(localStorage.getItem(DOCK_Y_KEY) ?? '')
    if (Number.isFinite(raw)) dockY.value = raw
  } catch { /* 忽略 */ }
  clampDockY()
  ro = new ResizeObserver(() => {
    pillW.value = pillEl.value?.offsetWidth ?? 0
    pillH.value = pillEl.value?.offsetHeight ?? 0
    clampDockY()
  })
  if (pillEl.value) ro.observe(pillEl.value)
  window.addEventListener('resize', clampDockY)
  // 欢迎词:挂载(入场幕布结束)即弹出,3s 后缩回缝隙,期间完全被动
  welcomeTimer = window.setTimeout(() => {
    welcome.value = false
    if (!bubbleJobId.value && ui.value === 'pop') ui.value = 'dock'
  }, 3000)
  timer = window.setInterval(() => {
    now.value = Date.now()
    if (ui.value === 'pop' && retractAt.value && Date.now() > retractAt.value) {
      ui.value = 'dock'
      bubbleJobId.value = null
    }
  }, 500)
})
onBeforeUnmount(() => {
  window.clearInterval(timer)
  window.clearTimeout(welcomeTimer)
  window.removeEventListener('resize', clampDockY)
  ro?.disconnect()
})

/* 未读事件 → 弹出一次(欢迎词让位给真实任务;面板开着时面板即提醒,不重复弹) */
watch(lastEvent, (evt) => {
  if (!evt || ui.value === 'panel') return
  welcome.value = false
  ui.value = 'pop'
  bubbleJobId.value = evt.jobId
  retractAt.value = Date.now() + RETRACT_MS
})

/* 面板翻开位置:上方放不下则翻到悬浮物下方 */
watch([ui, () => tasks.value.length, dockY], async () => {
  if (ui.value !== 'panel') return
  await nextTick()
  const h = panelEl.value?.offsetHeight ?? 0
  flipBelow.value = dockY.value - h - 44 < 0
})
</script>

<template>
  <div
    ref="launcherEl"
    class="tc-launcher"
    :class="[`st-${status}`, `md-${ui}`, { dragging }]"
    :style="{ top: dockY + 'px' }"
    @pointerdown="onDown"
    @pointermove="onMove"
    @pointerup="onUp"
    @pointerenter="onEnter"
    @pointerleave="onLeave"
  >
    <div ref="pillEl" class="tc-pill">
      <span class="tc-sliver">
        <span class="tc-bars"><i v-for="i in Math.min(unread, BAR_MAX)" :key="i"></i></span>
        <svg v-if="status === 'run'" :width="12" :height="Math.max(pillH, 34)"
             :viewBox="`0 0 12 ${Math.max(pillH, 34)}`" class="tc-svg">
          <path class="tc-march" fill="none" :style="{ stroke: 'var(--c-ink)' }"
                stroke-width="1.5" stroke-dasharray="25 25" pathLength="100" :d="sliverPath" />
        </svg>
      </span>
      <span class="tc-ring">
        <svg v-if="status === 'run'" :width="Math.max(pillW, 120)" :height="Math.max(pillH, 34)"
             :viewBox="`0 0 ${Math.max(pillW, 120)} ${Math.max(pillH, 34)}`" class="tc-svg">
          <path class="tc-march" fill="none" :style="{ stroke: 'var(--c-ink)' }"
                stroke-width="1.5" stroke-dasharray="25 25" pathLength="100" :d="ringPath" />
        </svg>
      </span>
      <span class="tc-pill-body">
        <svg v-if="status === 'run'" class="tc-ico" viewBox="0 0 18 18" width="16" height="16">
          <circle cx="9" cy="9" r="7" fill="none" stroke="var(--c-line)" stroke-width="2" />
          <circle class="tc-ico-arc" cx="9" cy="9" r="7" fill="none" stroke="var(--c-ink)"
                  stroke-width="2" stroke-linecap="round" stroke-dasharray="22 44" />
        </svg>
        {{ pillText }}
      </span>
    </div>

    <Transition name="tc-rise" appear
                appear-from-class="tc-rise-appear-from" appear-active-class="tc-rise-appear-active">
      <div v-if="welcome && !bubbleTask" class="tc-bubble">
        <div class="tc-bh">{{ t('web.jobs.welcome') }}</div>
      </div>
      <div v-else-if="bubbleTask" class="tc-bubble" :class="{ err: bubbleTask.status === 'error' }">
        <div class="tc-bh">
          {{ bubbleTask.status === 'error' ? '✗' : '✓' }}
          {{ bubbleTask.status === 'error'
            ? t('web.jobs.fail_long', { name: taskName(bubbleTask) })
            : t('web.jobs.done_long', { name: taskName(bubbleTask) }) }}
          <span class="tc-bx" @click="bubbleJobId = null">×</span>
        </div>
        <div class="tc-ba">
          <span class="tc-go" :class="{ err: bubbleTask.status === 'error' }" @click="openTask(bubbleTask)">
            {{ bubbleTask.status === 'error' ? t('web.jobs.view_log') : t('web.jobs.view_result') }}
          </span>
          <span v-if="bubbleTask.status === 'done'" class="tc-life">
            {{ t('web.jobs.keep_window', { t: fmt(remainingSec(bubbleTask)) }) }}
          </span>
        </div>
      </div>
    </Transition>

    <Transition name="tc-rise">
      <div v-if="ui === 'panel'" ref="panelEl" class="tc-panel" :class="{ below: flipBelow }">
        <div class="tc-ph">
          {{ t('web.jobs.title') }}
          <span class="tc-tog" @click="collapse">▴ {{ t('web.jobs.collapse') }}</span>
        </div>
        <div v-if="tasks.length === 0" class="tc-empty">{{ t('web.jobs.empty') }}</div>
        <div v-for="task in tasks" :key="task.jobId" class="tc-row">
          <template v-if="task.status === 'running'">
            <div class="tc-t1">
              {{ taskName(task) }}
              <span class="tc-st run">{{ t('web.jobs.running') }}</span>
            </div>
            <div class="tc-t2">
              <span class="tc-pbar"><i :style="{ width: (task.fraction ?? 0) * 100 + '%' }"></i></span>
              <span>{{ Math.round((task.fraction ?? 0) * 100) }}%</span>
            </div>
          </template>
          <template v-else-if="task.status === 'done'">
            <div class="tc-t1">
              {{ taskName(task) }}
              <span class="tc-st ok">✓ {{ t('web.jobs.done') }}</span>
            </div>
            <div class="tc-t2">
              <span class="tc-act" :class="{ off: remainingSec(task) <= 0 }" @click="openTask(task)">
                {{ t('web.jobs.view_result') }}
              </span>
              <span class="tc-cd">
                {{ remainingSec(task) > 0
                  ? t('web.jobs.keep_window', { t: fmt(remainingSec(task)) })
                  : t('web.jobs.archived') }}
              </span>
            </div>
          </template>
          <template v-else>
            <div class="tc-t1">
              {{ taskName(task) }}
              <span class="tc-st err">{{ t('web.jobs.fail') }}</span>
            </div>
            <div class="tc-t2">
              <span class="tc-errtext">{{ task.error === 'lost' ? t('web.jobs.lost') : task.error }}</span>
              <span class="tc-act" @click="openTask(task)">{{ t('web.jobs.view_log') }}</span>
            </div>
          </template>
        </div>
        <div class="tc-pf">{{ t('web.jobs.foot') }}</div>
      </div>
    </Transition>
  </div>
</template>

<style scoped>
/* ===== 悬浮物:贴右缘 fixed,位置(停靠高度)由 dockY 驱动 ===== */
.tc-launcher {
  position: fixed;
  right: 0;
  z-index: 80;
  cursor: grab;
  user-select: none;
  touch-action: none;
  /* 挂载入场:从右缘滑入(新插入元素上 CSS 动画自动播放,无需 Vue Transition) */
  animation: tc-in 0.7s cubic-bezier(0.22, 0.61, 0.36, 1);
}
@keyframes tc-in { from { opacity: 0; transform: translateX(26px); } }
@media (prefers-reduced-motion: reduce) { .tc-launcher { animation: none; } }
.tc-launcher.dragging { cursor: grabbing !important; }
.tc-pill {
  position: relative;
  display: inline-flex;
  align-items: center;
  padding: 8px 15px 8px 22px;
  font-size: 12.5px;
  white-space: nowrap;
  transition: transform 0.32s cubic-bezier(0.22, 0.61, 0.36, 1), box-shadow 0.2s ease;
  will-change: transform;
}
.tc-launcher.dragging .tc-pill { transition: none; }
.tc-pill-body {
  position: relative;
  z-index: 1;
  display: inline-flex;
  align-items: center;
  gap: 8px;
  transition: opacity 0.22s ease;
}
.tc-launcher.md-dock .tc-pill-body { opacity: 0; }
/* 缝隙:悬浮物左缘的"探头",停靠态可见 */
.tc-sliver {
  position: absolute;
  left: 0;
  top: 0;
  bottom: 0;
  width: 12px;
  border-radius: 8px 0 0 8px;
  background: var(--c-line);
  box-shadow: inset 0 0 0 1px var(--c-line-soft);
  overflow: hidden;
  transition: opacity 0.18s ease;
}
.tc-launcher.st-run .tc-sliver { background: var(--c-paper); }
.tc-launcher.st-unread .tc-sliver { background: var(--c-ink); }
.tc-launcher.md-pop .tc-sliver,
.tc-launcher.md-panel .tc-sliver { opacity: 0; }
.tc-launcher.md-dock:hover .tc-sliver { opacity: 0; }
/* 未读白条堆:1 条 = 1 个完成未读,封顶 5(黑底白条,任何底色上都可见) */
.tc-bars {
  position: absolute;
  inset: 0;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  gap: 2.5px;
}
.tc-bars i {
  width: 7px;
  height: 2.5px;
  border-radius: 1px;
  background: var(--c-paper);
}
/* 停靠:只露 12px;悬停带意图延迟后探头(扫过不探头,停住才探头) */
.tc-launcher.md-dock .tc-pill { transform: translateX(calc(100% - 12px)); }
.tc-launcher.md-dock:hover .tc-pill { transform: translateX(4px); transition-delay: 0.26s; }
.tc-launcher.md-dock:hover .tc-pill-body { opacity: 1; transition-delay: 0.26s; }
.tc-launcher.md-pop .tc-pill,
.tc-launcher.md-panel .tc-pill { transform: translateX(0); }
/* 三种身份对应的胶囊本体 */
.tc-launcher.st-run .tc-pill,
.tc-launcher.st-idle .tc-pill {
  background: var(--c-panel);
  border: 1px solid var(--c-line);
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.08);
  color: var(--c-mid);
}
.tc-launcher.st-unread .tc-pill {
  background: var(--c-ink);
  color: var(--c-on-accent);
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.14);
}
.tc-launcher.md-panel .tc-pill {
  background: var(--c-panel);
  color: var(--c-ink);
  border: 1px solid var(--c-line);
}
.tc-launcher.dragging .tc-pill { box-shadow: 0 6px 18px rgba(0, 0, 0, 0.22); }
/* 弹出/展开/悬停探头时运行中:轮转接到胶囊边框(SVG dashoffset,纯 CSS 动画) */
.tc-ring { display: none; position: absolute; inset: 0; pointer-events: none; }
.tc-svg { position: absolute; inset: 0; display: block; }
.tc-launcher.st-run.md-pop .tc-ring,
.tc-launcher.st-run.md-panel .tc-ring,
.tc-launcher.st-run.md-dock:hover .tc-ring { display: block; }
.tc-ico { display: block; flex: 0 0 auto; }
/* 运行指示圈:匀速旋转(转 = 有任务在跑;具体进度看面板各行与胶囊文字) */
.tc-ico-arc { transform-origin: 9px 9px; animation: tc-spin 1s linear infinite; }
@keyframes tc-spin { to { transform: rotate(360deg); } }
@media (prefers-reduced-motion: reduce) { .tc-ico-arc { animation: none; } }
/* 边框巡游:两黑两白(pathLength=100 归一化),偏移恰一个周期无缝循环 */
.tc-march { animation: tc-march 2.4s linear infinite; }
@keyframes tc-march { to { stroke-dashoffset: -50; } }
@media (prefers-reduced-motion: reduce) { .tc-march { animation: none; } }

/* ===== 锚定气泡(完成/失败弹出瞬间的反馈) ===== */
.tc-bubble {
  position: absolute;
  right: 10px;
  bottom: calc(100% + 12px);
  z-index: 2;
  width: 296px;
  background: var(--c-panel);
  border: 1px solid var(--c-line);
  border-left: 3px solid var(--c-ink);
  box-shadow: 0 4px 22px rgba(0, 0, 0, 0.12);
  padding: 12px 14px 11px;
}
.tc-bubble.err { border-left-color: var(--c-danger); }
.tc-bh {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 13px;
  font-weight: 600;
  overflow-wrap: anywhere;
}
.tc-bx { margin-left: auto; color: var(--c-ink3); cursor: pointer; font-weight: 400; font-size: 13px; padding-left: 10px; }
.tc-ba { display: flex; align-items: center; margin-top: 8px; gap: 12px; }
.tc-go {
  font-size: 12px;
  background: var(--c-ink);
  color: var(--c-on-accent);
  padding: 4px 13px;
  cursor: pointer;
  letter-spacing: 0.08em;
}
.tc-go.err { background: var(--c-danger); }
.tc-life { font-size: 11px; color: var(--c-ink3); }

/* ===== 展开面板(永远完整在界面内;行多时限高内部滚动) ===== */
.tc-panel {
  position: absolute;
  right: 10px;
  bottom: calc(100% + 12px);
  z-index: 2;
  width: 344px;
  max-height: min(70vh, 560px);
  overflow-y: auto;
  background: var(--c-panel);
  border: 1px solid var(--c-line);
  box-shadow: 0 4px 22px rgba(0, 0, 0, 0.1);
}
.tc-panel.below { top: calc(100% + 12px); bottom: auto; }
.tc-ph {
  display: flex;
  align-items: center;
  padding: 11px 16px;
  border-bottom: 1px solid var(--c-line-soft);
  font-size: 12px;
  letter-spacing: 0.2em;
  color: var(--c-mid);
  position: sticky;
  top: 0;
  background: var(--c-panel);
}
.tc-tog { margin-left: auto; color: var(--c-sub); cursor: pointer; font-size: 13px; letter-spacing: 0; }
.tc-row { padding: 11px 16px 12px; border-bottom: 1px solid var(--c-rowline); }
.tc-t1 { display: flex; align-items: center; gap: 8px; font-size: 13px; overflow-wrap: anywhere; }
.tc-st { margin-left: auto; font-size: 11px; letter-spacing: 0.08em; color: var(--c-sub); flex: 0 0 auto; }
.tc-st.run { color: var(--c-mid); }
.tc-st.ok { color: var(--c-ink); font-weight: 600; }
.tc-st.err { color: var(--c-danger); font-weight: 600; }
.tc-t2 { display: flex; align-items: center; margin-top: 8px; font-size: 11.5px; color: var(--c-dim); gap: 10px; }
.tc-pbar { flex: 1; height: 3px; background: var(--c-line-soft); position: relative; }
.tc-pbar i { position: absolute; inset: 0 auto 0 0; background: var(--c-ink); }
.tc-act {
  color: var(--c-ink);
  text-decoration: underline;
  text-underline-offset: 3px;
  cursor: pointer;
  flex: 0 0 auto;
}
.tc-act.off { color: var(--c-ink3); text-decoration: none; cursor: default; }
.tc-cd { flex: 0 0 auto; }
.tc-errtext { flex: 1; overflow-wrap: anywhere; }
.tc-pf {
  padding: 9px 16px;
  border-top: 1px solid var(--c-line-soft);
  font-size: 11px;
  color: var(--c-ink3);
  position: sticky;
  bottom: 0;
  background: var(--c-panel);
}
.tc-empty { padding: 18px 16px; font-size: 12px; color: var(--c-dim); }

/* ===== 入场/退场(rise 上浮淡入;180ms 淡出,避免硬消失) ===== */
.tc-rise-enter-active { transition: opacity 0.28s ease-out, transform 0.28s ease-out; }
.tc-rise-leave-active { transition: opacity 0.18s ease; }
.tc-rise-enter-from { opacity: 0; transform: translateY(6px); }
.tc-rise-leave-to { opacity: 0; }
/* 欢迎词首次入场:与胶囊本体同速(0.7s),交互动画(任务气泡/面板)保持 0.28s */
.tc-rise-appear-active { transition: opacity 0.7s cubic-bezier(0.22, 0.61, 0.36, 1), transform 0.7s cubic-bezier(0.22, 0.61, 0.36, 1); }
.tc-rise-appear-from { opacity: 0; transform: translateY(6px); }
</style>
