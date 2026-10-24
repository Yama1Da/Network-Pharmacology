<script setup lang="ts">
/** 运行日志(阶段 1):级别切换 + 按天日志文件清单与内容查看。
 * 骨架(顶栏/幕布/页脚)由 App 提供,样式见 styles/screen.css。 */
import { onMounted, ref } from 'vue'
import { api, type LogFile } from '../api'

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
    level.value = (await api.logsLevel()).level as 'INFO' | 'DEBUG'
  } catch {
    flash('日志级别读取失败')
  }
}

async function loadFiles(): Promise<void> {
  try {
    files.value = (await api.logsFiles()).files
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
    contents.value[name] = await api.logsContent(name)
  } catch {
    flash('日志内容读取失败')
  } finally {
    loadingContent.value = false
  }
}

async function setLevel(next: 'INFO' | 'DEBUG'): Promise<void> {
  if (next === level.value) return
  try {
    level.value = (await api.setLogsLevel(next)).level as 'INFO' | 'DEBUG'
    flash(`日志级别已切换为 ${level.value},对后续日志即时生效`)
  } catch {
    flash('级别切换失败')
  }
}

function refresh(): void {
  void loadFiles()
  void loadLevel()
}

const fmtSize = (n: number): string =>
  n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(1)} KB` : `${(n / 1048576).toFixed(1)} MB`

onMounted(() => {
  void loadLevel()
  void loadFiles()
})
</script>

<template>
  <p class="kicker s1">RUNTIME&nbsp;LOGS</p>
  <h1 class="title s2">运行日志</h1>
  <p class="sub s3">按天归档&nbsp;·&nbsp;保留两年&nbsp;·&nbsp;与分析进程共用同一份文件</p>

  <div class="panel s4">
    <div class="controls">
      <span class="ctl-label">级别</span>
      <span class="lv-group">
        <button class="lv" :class="{ on: level === 'INFO' }" @click="setLevel('INFO')">INFO</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: level === 'DEBUG' }" @click="setLevel('DEBUG')">DEBUG</button>
      </span>
      <i class="dot"></i>
      <span class="ctl-meta">当前 {{ level }} · 修改即时生效,重启后保留</span>
      <button class="go" title="刷新" @click="refresh">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <ul class="files">
      <li v-if="!files.length" class="empty">
        暂无日志文件 —— 跑一次分析后,这里会出现按天归档的记录。
      </li>
      <li v-for="f in files" :key="f.name" class="file">
        <button class="row as-btn" @click="toggle(f.name)">
          <span class="fname" :class="{ open: expanded === f.name }">
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
</template>

<style scoped>
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
.row.as-btn {
  width: 100%;
  border: none;
  background: none;
  font-family: inherit;
  cursor: pointer;
  text-align: left;
}
.fname {
  transition: color 0.25s ease;
}
.row:hover .fname {
  color: #000000;
}
.fname.open {
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
.content pre {
  margin: 0;
  padding: 10px 0 4px;
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
</style>
