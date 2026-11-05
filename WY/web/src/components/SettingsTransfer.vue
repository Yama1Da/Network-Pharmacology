<script setup lang="ts">
/** 数据迁移(设置页底部的导出/导入独立操作区,从 SettingsView 切出):
 * 导出走 job 打包后浏览器下载;导入两段式:上传解包预览 → 逐条决策 → 落盘。
 * 密钥绝不进包:导出内容只有结果目录与 ui_settings.json(白名单设置,无密钥)。
 * 与宿主的接口:flash 顶栏提示、reload 全量刷新快照、ui-settings 导入后经
 * ui-imported 事件交宿主回填主题/字号/语言(宿主的 applyImportedUi)。 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot, type TransferImportReport, type TransferItemReport } from '../api'
import { useJob } from '../composables/useJob'
import { downloadBlob } from '../utils/download'

const props = defineProps<{
  /** 顶部 notice 的提示函数(宿主 useNotice 的 flash) */
  flash: (msg: string) => void
  /** 全量刷新设置快照(宿主的 load) */
  reload: () => Promise<void>
}>()

const emit = defineEmits<{ (e: 'ui-imported', snap: SettingsSnapshot): void }>()

const { t } = useI18n()



type ScopeKey = 'analysis' | 'star' | 'geo' | 'ui_settings'
const SCOPE_DEFS: { key: ScopeKey; label: string }[] = [
  { key: 'analysis', label: 'web.transfer.scope_analysis' },
  { key: 'star', label: 'web.transfer.scope_star' },
  { key: 'geo', label: 'web.transfer.scope_geo' },
  { key: 'ui_settings', label: 'web.transfer.scope_settings' }
]
const KIND_LABEL: Record<string, string> = Object.fromEntries(
  SCOPE_DEFS.map((s) => [s.key, s.label])
)
const exportScopes = ref<Record<ScopeKey, boolean>>({
  analysis: true, star: true, geo: true, ui_settings: true
})

const exportJob = useJob((result) => {
  const filename = String(result.filename ?? '')
  void (async () => {
    try {
      downloadBlob(await api.transferDownload(filename), filename)
      const counts = (result.counts ?? {}) as Record<string, number>
      const total = Object.values(counts).reduce((a, b) => a + Number(b), 0)
      props.flash(t('web.transfer.export_done', { n: total }))
    } catch {
      props.flash(t('web.transfer.download_fail'))
    }
  })()
})

async function runExport(): Promise<void> {
  const scopes = SCOPE_DEFS.filter((s) => exportScopes.value[s.key]).map((s) => s.key)
  if (!scopes.length) {
    props.flash(t('web.transfer.no_scope'))
    return
  }
  await exportJob.run(api.transferExport(scopes), { failText: t('web.transfer.export_fail') })
}

/* 导入两段式:预览行 → 决策表(action + 改名输入) → commit;失败项留在会话可重试 */
interface RowDecision {
  action: 'import' | 'overwrite' | 'skip' | 'rename'
  new_name: string
}
const importing = ref(false)
const importSid = ref('')
const importReport = ref<TransferImportReport | null>(null)
const decisions = ref<Record<string, RowDecision>>({})
const uiDone = ref(false)

const keyOf = (kind: string, name: string): string => `${kind}/${name}`
const selectedCount = computed(
  () => Object.values(decisions.value).filter((d) => d.action !== 'skip').length
)

function dec(kind: string, name: string): RowDecision {
  const key = keyOf(kind, name)
  let d = decisions.value[key]
  if (!d) {
    d = { action: 'skip', new_name: '' }
    decisions.value[key] = d
  }
  return d
}

function filePicked(e: Event): void {
  const input = e.target as HTMLInputElement
  const f = input.files?.[0]
  input.value = '' // 允许重复选择同一个文件重新解析
  if (f) void uploadImport(f)
}

async function uploadImport(f: File): Promise<void> {
  importing.value = true
  resetImport()
  try {
    const r = await api.transferImport(f)
    importSid.value = r.sid
    importReport.value = r.report
    // 默认决策:无冲突 → 导入;与现有数据重名或不可导入 → 跳过(由用户显式改选)
    for (const it of r.report.items) {
      dec(it.kind, it.name).action = it.conflict || it.invalid ? 'skip' : 'import'
    }
    if (r.report.ui_settings.present) {
      // 界面设置无害(白名单过滤,覆盖的只是阈值/主题类偏好),默认直接接入,
      // 重名只作提示;结果目录的冲突才默认跳过(覆盖会动真实数据)
      dec('ui_settings', 'ui_settings.json').action = 'import'
    }
  } catch (e) {
    props.flash(e instanceof Error ? e.message : t('web.transfer.import_fail'))
  } finally {
    importing.value = false
  }
}

async function commitImport(): Promise<void> {
  if (!importSid.value) return
  importing.value = true
  try {
    // skip 决策也发后端:由 skipped 回执统一收账,本地不做过滤
    const payload = Object.entries(decisions.value).map(([key, d]) => {
      const idx = key.indexOf('/')
      return {
        kind: key.slice(0, idx),
        name: key.slice(idx + 1),
        action: d.action === 'import' ? ('overwrite' as const) : d.action,
        new_name: d.new_name
      }
    })
    const r = await api.transferCommit(importSid.value, payload)
    const rep = r.report
    const parts = [
      rep.applied.length ? t('web.transfer.done_applied', { n: rep.applied.length }) : '',
      rep.skipped.length ? t('web.transfer.done_skipped', { n: rep.skipped.length }) : '',
      rep.failed.length ? t('web.transfer.done_failed', { n: rep.failed.length }) : ''
    ].filter(Boolean)
    props.flash(parts.join(' · ') || t('web.transfer.import_none'))
    if (r.ui_snapshot) emit('ui-imported', r.ui_snapshot)
    // 已落定(导入/跳过)的行从预览里移除,失败行留下可改决策重试
    const done = new Set(
      [...rep.applied, ...rep.skipped].map((a) => keyOf(a.kind, a.name))
    )
    if (importReport.value) {
      importReport.value = {
        ...importReport.value,
        items: importReport.value.items.filter((it) => !done.has(keyOf(it.kind, it.name)))
      }
      if (rep.applied.some((a) => a.kind === 'ui_settings')) uiDone.value = true
    }
    decisions.value = Object.fromEntries(
      Object.entries(decisions.value).filter(([k]) => !done.has(k))
    )
    if (!rep.failed.length) resetImport()
    await props.reload()
  } catch (e) {
    props.flash(e instanceof Error ? e.message : t('web.transfer.import_fail'))
  } finally {
    importing.value = false
  }
}

async function cancelImport(): Promise<void> {
  const sid = importSid.value
  resetImport()
  if (sid) {
    try {
      await api.transferCancel(sid)
    } catch {
      /* 会话可能已过期,忽略 */
    }
  }
}

function resetImport(): void {
  importSid.value = ''
  importReport.value = null
  decisions.value = {}
  uiDone.value = false
}

function humanSize(n: number): string {
  if (n >= 1 << 30) return `${(n / (1 << 30)).toFixed(1)} GB`
  if (n >= 1 << 20) return `${(n / (1 << 20)).toFixed(1)} MB`
  if (n >= 1 << 10) return `${(n / (1 << 10)).toFixed(1)} KB`
  return `${n} B`
}

const kindLabel = (it: TransferItemReport): string => t(KIND_LABEL[it.kind] ?? it.kind)
</script>


<template>
  <!-- ======== 数据迁移:导出 / 两段式导入(独立操作区,与上方保存按钮无关) ======== -->
  <section class="sec">
    <h3 class="sec-title">{{ t('web.transfer.title') }}</h3>
    <p class="hint-line">{{ t('web.transfer.hint') }}</p>

    <div class="tf-scopes">
      <label v-for="s in SCOPE_DEFS" :key="s.key" class="chk">
        <input v-model="exportScopes[s.key]" type="checkbox" /> {{ t(s.label) }}
      </label>
    </div>
    <button class="bar-btn" :disabled="exportJob.running.value || importing"
            @click="runExport">{{ t('web.transfer.export_btn') }}</button>
    <div v-if="exportJob.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((exportJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ exportJob.text.value || t('web.transfer.exporting') }}</p>
    </div>
    <p v-if="exportJob.error.value" class="notice">{{ exportJob.error.value }}</p>

    <div class="tf-import-row">
      <label class="bar-btn tf-file-btn">
        <input type="file" accept=".zip" @change="filePicked" />
        {{ t('web.transfer.import_btn') }}
      </label>
      <span v-if="importing && !importReport" class="hint-line">
        {{ t('web.transfer.import_analyzing') }}
      </span>
    </div>

    <template v-if="importReport">
      <p v-for="w in importReport.warnings" :key="w" class="notice">{{ w }}</p>
      <div class="tf-table">
        <div class="tf-row tf-head">
          <span>{{ t('web.transfer.col_kind') }}</span>
          <span>{{ t('web.transfer.col_name') }}</span>
          <span>{{ t('web.transfer.col_size') }}</span>
          <span>{{ t('web.transfer.col_action') }}</span>
        </div>
        <div v-for="it in importReport.items" :key="keyOf(it.kind, it.name)" class="tf-row">
          <span class="tf-kind">{{ kindLabel(it) }}</span>
          <span class="tf-name">
            {{ it.name }}
            <em v-if="it.conflict" class="tf-flag">{{ t('web.transfer.conflict') }}</em>
            <em v-if="it.invalid" class="tf-flag bad" :title="it.invalid_reason">
              {{ t('web.transfer.invalid_item') }}
            </em>
          </span>
          <span class="tf-size">{{ humanSize(it.bytes) }}</span>
          <span class="tf-action">
            <select v-model="dec(it.kind, it.name).action" class="tf-select"
                    :disabled="it.invalid">
              <option value="import">{{ t('web.transfer.act_import') }}</option>
              <option v-if="it.conflict" value="overwrite">{{ t('web.transfer.act_overwrite') }}</option>
              <option value="skip">{{ t('web.transfer.act_skip') }}</option>
              <option v-if="it.conflict && !it.invalid" value="rename">{{ t('web.transfer.act_rename') }}</option>
            </select>
            <input v-if="dec(it.kind, it.name).action === 'rename'"
                   v-model="dec(it.kind, it.name).new_name" class="tf-new"
                   :placeholder="t('web.transfer.new_name_ph')" spellcheck="false" />
          </span>
        </div>
        <div v-if="importReport.ui_settings.present && !uiDone" class="tf-row">
          <span class="tf-kind">{{ t('web.transfer.scope_settings') }}</span>
          <span class="tf-name">ui_settings.json
            <em v-if="importReport.ui_settings.conflict" class="tf-flag">
              {{ t('web.transfer.conflict') }}
            </em>
          </span>
          <span class="tf-size">—</span>
          <span class="tf-action">
            <select v-model="dec('ui_settings', 'ui_settings.json').action" class="tf-select">
              <option value="import">{{ t('web.transfer.act_import') }}</option>
              <option value="skip">{{ t('web.transfer.act_skip') }}</option>
            </select>
          </span>
        </div>
      </div>
      <div class="ver-row">
        <button class="bar-btn" :disabled="importing || !selectedCount"
                @click="commitImport">{{ t('web.transfer.run_import') }}</button>
        <button class="bar-btn tf-ghost" :disabled="importing"
                @click="cancelImport">{{ t('web.transfer.cancel_import') }}</button>
      </div>
    </template>
  </section>
</template>

<style scoped>
/* 以下 .sec/.sec-title/.bar-btn/.ver-row/.notice 为宿主共享类的本地副本
(scoped 样式不跨组件,子组件需自持;若日后拆出更多设置子组件,应上提共享) */
.sec {
  margin-top: 26px;
  padding-top: 18px;
  border-top: 1px solid var(--c-line-soft);
}
.sec-title {
  margin: 0 0 8px;
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0.14em;
  color: var(--c-ink);
}
.bar-btn {
  margin-top: 18px;
  border: none;
  background: none;
  font-family: inherit;
  font-size: 12px;
  letter-spacing: 0.08em;
  color: var(--c-ink);
  cursor: pointer;
  padding: 4px 0;
  border-bottom: 1px solid var(--c-accent);
  transition: opacity 0.25s ease;
}
.bar-btn:disabled {
  cursor: default;
  border-bottom-color: var(--c-ghost);
  color: var(--c-ink4);
}
.ver-row {
  display: flex;
  align-items: baseline;
  gap: 12px;
  padding: 6px 0;
  flex-wrap: wrap;
}
.notice {
  margin: 0 0 14px;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: var(--c-dim);
}

/* ======== 数据迁移(导出/导入) ======== */
.tf-scopes {
  display: flex;
  flex-wrap: wrap;
  gap: 6px 28px;
  margin-top: 10px;
}
.tf-import-row {
  display: flex;
  align-items: baseline;
  gap: 18px;
  margin-top: 22px;
}
.tf-file-btn {
  margin-top: 0;
}
.tf-file-btn input[type="file"] {
  display: none;
}
.tf-table {
  margin-top: 14px;
  border-top: 1px solid var(--c-line-soft);
}
.tf-row {
  display: grid;
  grid-template-columns: 92px minmax(0, 1fr) 74px minmax(180px, auto);
  gap: 12px;
  align-items: center;
  padding: 8px 0;
  border-bottom: 1px solid var(--c-line-soft);
  font-size: 12.5px;
}
.tf-row.tf-head {
  font-size: 10.5px;
  letter-spacing: 0.14em;
  color: var(--c-ink3);
}
.tf-kind {
  color: var(--c-ink2);
}
.tf-name {
  min-width: 0;
  overflow-wrap: anywhere;
}
.tf-flag {
  display: inline-block;
  margin-left: 8px;
  padding: 1px 6px;
  border: 1px solid var(--c-line);
  font-style: normal;
  font-size: 10.5px;
  letter-spacing: 0.06em;
  color: var(--c-ink2);
}
.tf-flag.bad {
  color: var(--c-danger);
  border-color: var(--c-danger);
}
.tf-size {
  font-variant-numeric: tabular-nums;
  color: var(--c-ink3);
  text-align: right;
}
.tf-action {
  display: flex;
  gap: 8px;
  align-items: center;
}
.tf-select,
.tf-new {
  border: none;
  border-bottom: 1px solid var(--c-line);
  background: transparent;
  font-family: inherit;
  font-size: 12.5px;
  color: var(--c-ink);
  padding: 3px 2px;
  outline: none;
}
.tf-new {
  width: 130px;
}
.tf-ghost {
  border-bottom-color: transparent;
  color: var(--c-ink4);
}
.tf-ghost:hover {
  color: var(--c-ink);
}
.tf-ghost:disabled {
  color: var(--c-ink4);
}
</style>
