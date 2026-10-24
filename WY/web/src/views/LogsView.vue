<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { api, type LogFile } from '../api'

/* ---------------- 数据 ---------------- */

const level = ref<'INFO' | 'DEBUG'>('INFO')
const files = ref<LogFile[]>([])
const expanded = ref('')
const contents = ref<Record<string, { content: string; truncated: boolean }>>({})
const loadingContent = ref(false)
const notice = ref('')
let noticeTimer = 0

function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 2600)
}

async function loadLevel(): Promise<void> {
  try {
    level.value = (await api.level()).level as 'INFO' | 'DEBUG'
  } catch {
    flash('日志级别读取失败')
  }
}

async function loadFiles(): Promise<void> {
  try {
    files.value = (await api.files()).files
  } catch {
    flash('日志列表读取失败')
  }
}

async function toggle(name: string): Promise<void> {
  expanded.value = expanded.value === name ? '' : name
  if (expanded.value && !contents.value[name]) await loadContent(name)
}

async function loadContent(name: string): Promise<void> {
  loadingContent.value = true
  try {
    contents.value[name] = await api.content(name)
  } catch {
    flash('日志内容读取失败')
  } finally {
    loadingContent.value = false
  }
}

async function setLevel(next: 'INFO' | 'DEBUG'): Promise<void> {
  if (next === level.value) return
  try {
    level.value = (await api.setLevel(next)).level as 'INFO' | 'DEBUG'
    flash(`日志级别已切换为 ${level.value},对后续日志即时生效`)
  } catch {
    flash('级别切换失败')
  }
}

function refresh(): void {
  void loadFiles()
  void loadLevel()
}

const totalSize = computed(() => files.value.reduce((s, f) => s + f.size, 0))
const fmtSize = (n: number): string =>
  n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(1)} KB` : `${(n / 1048576).toFixed(1)} MB`

/* ---------------- 入场动画状态(与主页同款) ---------------- */

const introVisible = ref(true)
const ready = ref(false)
let introTimer = 0

onMounted(() => {
  void loadLevel()
  void loadFiles()
  introTimer = window.setTimeout(() => {
    introVisible.value = false
    ready.value = true
  }, 1750)
})

onBeforeUnmount(() => {
  window.clearTimeout(noticeTimer)
  window.clearTimeout(introTimer)
})
</script>

<template>
  <div class="screen-logs">
    <!-- 背景:两圈几乎不可见的呼吸细环 -->
    <div class="whisper w2" aria-hidden="true"></div>
    <div class="whisper w1" aria-hidden="true"></div>

    <!-- 入场幕布:细环自绘 + 字标,随后淡出揭幕 -->
    <Transition name="intro-fade">
      <div v-if="introVisible" class="intro" aria-hidden="true">
        <svg class="intro-ring" viewBox="0 0 64 64" width="56" height="56">
          <circle class="intro-ring-path" cx="32" cy="32" r="30" />
        </svg>
        <span class="intro-word">NPAP</span>
      </div>
    </Transition>

    <!-- 顶栏 -->
    <header class="top">
      <span class="brand">NPAP<i class="brand-dot"></i></span>
      <nav class="links">
        <a href="http://127.0.0.1:8501" title="过渡期:打开 Streamlit 主界面">主界面</a>
        <a @click.prevent="refresh">刷新</a>
      </nav>
    </header>

    <!-- 中央 -->
    <main class="center" :class="{ ready }">
      <p class="kicker s1">RUNTIME&nbsp;LOGS</p>
      <h1 class="title s2">运行日志</h1>
      <p class="sub s3">按天归档&nbsp;·&nbsp;保留两年&nbsp;·&nbsp;与分析进程共用同一份文件</p>

      <div class="panel s4">
        <!-- 级别切换:文本开关,激活者着墨色 -->
        <div class="controls">
          <span class="ctl-label">级别</span>
          <span class="lv-group">
            <button class="lv" :class="{ on: level === 'INFO' }" @click="setLevel('INFO')">INFO</button>
            <i class="lv-sep">/</i>
            <button class="lv" :class="{ on: level === 'DEBUG' }" @click="setLevel('DEBUG')">DEBUG</button>
          </span>
          <i class="dot"></i>
          <span class="ctl-meta">当前 {{ level }} · 修改即时生效,重启后保留</span>
          <button class="go" title="刷新" aria-label="刷新" @click="refresh">
            <svg viewBox="0 0 20 20" width="16" height="16">
              <path
                d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4"
                fill="none"
                stroke="currentColor"
                stroke-width="1.6"
                stroke-linecap="round"
                stroke-linejoin="round"
              />
            </svg>
          </button>
        </div>

        <Transition name="note">
          <p v-if="notice" class="notice">{{ notice }}</p>
        </Transition>

        <!-- 文件列表:发丝线分隔,点击展开内容 -->
        <ul class="files">
          <li v-if="!files.length" class="empty">
            暂无日志文件 —— 跑一次分析后,这里会出现按天归档的记录。
          </li>
          <li
            v-for="f in files"
            :key="f.name"
            class="file"
            :class="{ open: expanded === f.name }"
          >
            <button class="row" @click="toggle(f.name)">
              <span class="fname">
                <i v-if="f.is_current" class="cur-dot" title="当前正在写入"></i>{{ f.name }}
              </span>
              <span class="fmeta">{{ fmtSize(f.size) }} · {{ f.modified }}</span>
            </button>
            <div v-if="expanded === f.name" class="content">
              <p v-if="loadingContent && !contents[f.name]" class="hint-line">读取中…</p>
              <pre v-else>{{ contents[f.name]?.content }}</pre>
              <p v-if="contents[f.name]?.truncated" class="trunc">
                文件过大,仅显示末尾 2 MB —— 完整内容请直接打开日志文件。
              </p>
            </div>
          </li>
        </ul>
      </div>
    </main>

    <!-- 底部 -->
    <footer class="foot">
      <span class="foot-stats">
        <b>{{ files.length }}</b> 个文件 · 共 <b>{{ fmtSize(totalSize) }}</b>
      </span>
      <i class="foot-dot"></i>
      <span>© 2026 NPAP · 仅供科研使用</span>
    </footer>
  </div>
</template>

<style scoped>
.screen-logs {
  position: relative;
  width: 100%;
  height: 100%;
  overflow: hidden;
  background: #fdfdfd;
  color: #1a1a1a;
  font-family: 'Helvetica Neue', Helvetica, Arial, 'Microsoft YaHei', sans-serif;
  display: flex;
  flex-direction: column;
}

/* ---------- 入场幕布(与主页同款) ---------- */
.intro {
  position: absolute;
  inset: 0;
  z-index: 50;
  background: #fdfdfd;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 22px;
}
.intro-ring-path {
  fill: none;
  stroke: #1a1a1a;
  stroke-width: 1;
  stroke-dasharray: 188.5;
  animation: ring-draw 1.15s cubic-bezier(0.4, 0, 0.2, 1) 0.15s both;
}
@keyframes ring-draw {
  from {
    stroke-dashoffset: 188.5;
  }
  to {
    stroke-dashoffset: 0;
  }
}
.intro-word {
  font-size: 11px;
  font-weight: 600;
  letter-spacing: 0.55em;
  text-indent: 0.55em;
  color: #1a1a1a;
  animation: word-in 1s ease 0.35s both;
}
@keyframes word-in {
  from {
    opacity: 0;
    transform: translateY(4px);
  }
  to {
    opacity: 1;
    transform: translateY(0);
  }
}
.intro-fade-leave-active {
  transition: opacity 0.65s ease;
}
.intro-fade-leave-to {
  opacity: 0;
}

/* ---------- 呼吸细环 ---------- */
.whisper {
  position: absolute;
  top: 50%;
  left: 50%;
  border: 1px solid #f0f0ee;
  border-radius: 50%;
  transform: translate(-50%, -54%);
  pointer-events: none;
}
.w1 {
  width: min(76vmin, 720px);
  height: min(76vmin, 720px);
  animation: breathe 11s ease-in-out infinite alternate;
}
.w2 {
  width: min(88vmin, 860px);
  height: min(88vmin, 860px);
  border-color: #f6f6f4;
  animation: breathe 15s ease-in-out infinite alternate-reverse;
}
@keyframes breathe {
  to {
    transform: translate(-50%, -54%) scale(1.035);
  }
}

/* ---------- 顶栏 ---------- */
.top {
  position: relative;
  z-index: 2;
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 28px 48px;
}
.brand {
  font-size: 13px;
  font-weight: 600;
  letter-spacing: 0.3em;
}
.brand-dot {
  display: inline-block;
  width: 4px;
  height: 4px;
  background: #1a1a1a;
  margin-left: 7px;
  vertical-align: 2px;
}
.links {
  display: flex;
  gap: 30px;
}
.links a {
  font-size: 11.5px;
  letter-spacing: 0.14em;
  color: #a6a6a6;
  cursor: pointer;
  transition: color 0.25s ease;
  user-select: none;
}
.links a:hover {
  color: #1a1a1a;
}

/* ---------- 中央 ---------- */
.center {
  position: relative;
  z-index: 2;
  flex: 1;
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 34px 24px 40px;
}

.s1,
.s2,
.s3,
.s4,
.s5 {
  opacity: 0;
}
.ready .s1,
.ready .s2,
.ready .s3,
.ready .s4,
.ready .s5 {
  animation: fadeUp 0.8s cubic-bezier(0.22, 0.61, 0.36, 1) both;
}
.ready .s2 {
  animation-delay: 0.08s;
}
.ready .s3 {
  animation-delay: 0.16s;
}
.ready .s4 {
  animation-delay: 0.26s;
}
.ready .s5 {
  animation-delay: 0.36s;
}
@keyframes fadeUp {
  from {
    opacity: 0;
    transform: translateY(10px);
  }
  to {
    opacity: 1;
    transform: translateY(0);
  }
}

.kicker {
  margin: 0 0 22px;
  font-size: 10px;
  letter-spacing: 0.62em;
  text-indent: 0.62em;
  color: #c2c2c2;
  text-align: center;
}
.title {
  margin: 0;
  font-size: clamp(28px, 3.4vw, 40px);
  font-weight: 300;
  letter-spacing: 0.18em;
  text-indent: 0.18em;
  color: #1a1a1a;
  text-align: center;
}
.sub {
  margin: 16px 0 0;
  font-size: 12px;
  letter-spacing: 0.34em;
  text-indent: 0.34em;
  color: #adadad;
  text-align: center;
}

/* ---------- 级别控制行:发丝线基线,控件坐在线上 ---------- */
.panel {
  width: min(760px, 92vw);
  margin-top: 48px;
}
.controls {
  display: flex;
  align-items: baseline;
  gap: 16px;
  padding-bottom: 12px;
  border-bottom: 1px solid #e6e6e4;
}
.ctl-label {
  font-size: 11px;
  letter-spacing: 0.3em;
  color: #c2c2c2;
}
.lv-group {
  display: inline-flex;
  align-items: baseline;
  gap: 10px;
}
.lv {
  border: none;
  background: none;
  padding: 2px 0;
  font-family: inherit;
  font-size: 12.5px;
  letter-spacing: 0.12em;
  color: #cfcfcd;
  cursor: pointer;
  transition: color 0.25s ease;
}
.lv:hover {
  color: #8a8a88;
}
.lv.on {
  color: #1a1a1a;
  font-weight: 500;
}
.lv-sep {
  color: #dedede;
  font-size: 12px;
}
.dot {
  width: 3px;
  height: 3px;
  border-radius: 50%;
  background: #dedede;
  align-self: center;
}
.ctl-meta {
  font-size: 11px;
  letter-spacing: 0.06em;
  color: #bcbcbc;
}
.go {
  margin-left: auto;
  align-self: center;
  width: 38px;
  height: 38px;
  border: none;
  border-radius: 50%;
  background: #1a1a1a;
  color: #fdfdfd;
  display: grid;
  place-items: center;
  cursor: pointer;
  font-family: inherit;
  transition: transform 0.3s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.go:hover {
  transform: scale(1.08);
}
.go:active {
  transform: scale(0.94);
}

/* ---------- 无框文字提示 ---------- */
.notice {
  margin: 14px 0 0;
  font-size: 11px;
  letter-spacing: 0.1em;
  color: #9a9a98;
}
.note-enter-active,
.note-leave-active {
  transition: all 0.4s ease;
}
.note-enter-from,
.note-leave-to {
  opacity: 0;
  transform: translateY(4px);
}

/* ---------- 文件列表:发丝线分隔行 ---------- */
.files {
  list-style: none;
  margin: 0;
  padding: 0;
}
.file {
  border-bottom: 1px solid #f0f0ee;
}
.row {
  width: 100%;
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 16px;
  border: none;
  background: none;
  padding: 14px 2px;
  font-family: inherit;
  cursor: pointer;
  text-align: left;
}
.fname {
  font-size: 13px;
  letter-spacing: 0.04em;
  color: #6e6e6c;
  font-variant-numeric: tabular-nums;
  transition: color 0.25s ease;
}
.row:hover .fname {
  color: #1a1a1a;
}
.file.open .fname {
  color: #1a1a1a;
  font-weight: 500;
}
.cur-dot {
  display: inline-block;
  width: 4px;
  height: 4px;
  border-radius: 50%;
  background: #1a1a1a;
  margin-right: 9px;
  vertical-align: 2.5px;
}
.fmeta {
  font-size: 11px;
  letter-spacing: 0.08em;
  color: #bcbcbc;
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}
.content {
  padding: 0 2px 18px;
}
.content pre {
  margin: 0;
  padding: 12px 0 4px;
  max-height: 420px;
  overflow: auto;
  font-family: 'SF Mono', 'Cascadia Code', Consolas, 'Courier New', monospace;
  font-size: 11.5px;
  line-height: 1.75;
  color: #6e6e6c;
  white-space: pre-wrap;
  word-break: break-all;
  tab-size: 4;
}
.hint-line {
  margin: 10px 0 0;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: #c6c6c4;
}
.trunc {
  margin: 10px 0 0;
  font-size: 10.5px;
  letter-spacing: 0.06em;
  color: #c6c6c4;
}
.empty {
  padding: 30px 2px;
  font-size: 12px;
  letter-spacing: 0.06em;
  color: #c6c6c4;
  text-align: center;
}

/* ---------- 底部 ---------- */
.foot {
  position: relative;
  z-index: 2;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 16px;
  padding: 20px 24px;
  border-top: 1px solid #f0f0ee;
  font-size: 10.5px;
  letter-spacing: 0.1em;
  color: #bcbcbc;
}
.foot-stats b {
  font-weight: 500;
  color: #a3a3a3;
  font-variant-numeric: tabular-nums;
}
.foot-dot {
  width: 3px;
  height: 3px;
  border-radius: 50%;
  background: #dedede;
}

/* ---------- 弱动效偏好 ---------- */
@media (prefers-reduced-motion: reduce) {
  .intro {
    display: none;
  }
  .s1,
  .s2,
  .s3,
  .s4,
  .s5,
  .ready .s1,
  .ready .s2,
  .ready .s3,
  .ready .s4,
  .ready .s5 {
    animation: none;
    opacity: 1;
  }
  .w1,
  .w2 {
    animation: none;
  }
}

@media (max-width: 760px) {
  .top {
    padding: 20px 26px;
  }
  .controls {
    flex-wrap: wrap;
  }
  .row {
    flex-direction: column;
    gap: 4px;
  }
}
</style>
