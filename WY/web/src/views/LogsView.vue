<script setup lang="ts">
/** 运行日志(阶段 1):级别切换 + 按天日志文件清单与内容查看。
 * 骨架(顶栏/幕布/页脚)由 App 提供,样式见 styles/screen.css。 */
import { onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type LogFile } from '../api'
import { useNotice } from '../composables/useNotice'
import { fmtSize } from '../utils/format'

const { t } = useI18n()
const { notice, flash } = useNotice(2600)
const level = ref<'INFO' | 'DEBUG'>('INFO')
const files = ref<LogFile[]>([])
const expanded = ref('')
const contents = ref<Record<string, { content: string; truncated: boolean }>>({})
const loadingContent = ref(false)

async function loadLevel(): Promise<void> {
  try {
    level.value = (await api.logsLevel()).level as 'INFO' | 'DEBUG'
  } catch {
    flash(t('web.logs.level_read_fail'))
  }
}

async function loadFiles(): Promise<void> {
  try {
    files.value = (await api.logsFiles()).files
  } catch {
    flash(t('web.logs.list_read_fail'))
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
    flash(t('web.logs.content_read_fail'))
  } finally {
    loadingContent.value = false
  }
}

async function setLevel(next: 'INFO' | 'DEBUG'): Promise<void> {
  if (next === level.value) return
  try {
    level.value = (await api.setLogsLevel(next)).level as 'INFO' | 'DEBUG'
    flash(t('web.logs.level_set', { level: level.value }))
  } catch {
    flash(t('web.logs.level_set_fail'))
  }
}

function refresh(): void {
  void loadFiles()
  void loadLevel()
}

onMounted(() => {
  void loadLevel()
  void loadFiles()
})
</script>

<template>
  <p class="kicker s1">{{ t('web.logs.kicker') }}</p>
  <h1 class="title s2">{{ t('web.logs.title') }}</h1>
  <p class="sub s3">{{ t('web.logs.sub') }}</p>

  <div class="panel s4">
    <div class="controls">
      <span class="ctl-label">{{ t('web.logs.level_label') }}</span>
      <span class="lv-group">
        <button class="lv" :class="{ on: level === 'INFO' }" @click="setLevel('INFO')">INFO</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: level === 'DEBUG' }" @click="setLevel('DEBUG')">DEBUG</button>
      </span>
      <i class="dot"></i>
      <span class="ctl-meta">{{ t('web.logs.level_meta', { level }) }}</span>
      <button class="go" :title="t('web.common.refresh')" @click="refresh">
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
      <li v-if="!files.length" class="empty">{{ t('web.logs.empty') }}</li>
      <li v-for="f in files" :key="f.name" class="file">
        <button class="row as-btn" @click="toggle(f.name)">
          <span class="fname" :class="{ open: expanded === f.name }">
            <i v-if="f.is_current" class="cur-dot" :title="t('web.logs.current_dot')"></i>{{ f.name }}
          </span>
          <span class="fmeta">{{ fmtSize(f.size) }} · {{ f.modified }}</span>
        </button>
        <div v-if="expanded === f.name" class="content">
          <p v-if="loadingContent && !contents[f.name]" class="hint-line">{{ t('web.logs.reading') }}</p>
          <pre v-else>{{ contents[f.name]?.content }}</pre>
          <p v-if="contents[f.name]?.truncated" class="trunc">{{ t('web.logs.truncated') }}</p>
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
  color: var(--c-faint);
  cursor: pointer;
  transition: color 0.25s ease;
}
.lv:hover {
  color: var(--c-ink2);
}
.lv.on {
  color: var(--c-ink);
  font-weight: 500;
}
.lv-sep {
  color: var(--c-ghost);
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
  color: var(--c-ink);
}
.fname.open {
  font-weight: 500;
}
.cur-dot {
  display: inline-block;
  width: 4px;
  height: 4px;
  border-radius: 50%;
  background: var(--c-accent);
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
  color: var(--c-mid);
  white-space: pre-wrap;
  word-break: break-all;
  tab-size: 4;
}
.hint-line {
  margin: 10px 0 0;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: var(--c-ink4);
}
.trunc {
  margin: 10px 0 0;
  font-size: 10.5px;
  letter-spacing: 0.06em;
  color: var(--c-ink4);
}
</style>
