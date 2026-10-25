<script setup lang="ts">
/** ML 结果管理(阶段 3c,ML_Results_Information 的 TS 迁移,独立页面):
 * 全部 ML run 清单(双根同扫:GEOresult + results/mlscreen)· 关键词筛选 ·
 * 勾选删除(两段确认)· ZIP 下载 · keep-N 清理。数据层在 core.ml_results。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from '../api'
import Stepper from '../components/Stepper.vue'

const { t } = useI18n()
const runs = ref<Record<string, any>[]>([])
const runsLoaded = ref(false)
const runsQuery = ref('')
const runsSel = ref<Set<string>>(new Set())
const runsPage = ref(1)
const RUNS_PAGE = 10
const cleanupKeep = ref(6)
const cleanupBusy = ref(false)
const cleanupErr = ref('')
const runsDeleteBusy = ref(false)
const runsDeleteErr = ref('')
const runsConfirm = ref(false)
let runsConfirmTimer = 0
const zipBusy = ref(false)

const notice = ref('')
let noticeTimer = 0
function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 4200)
}

async function loadRuns(): Promise<void> {
  try {
    runs.value = (await api.mlRuns()).runs
    runsLoaded.value = true
    const names = new Set(runs.value.map((r) => r.name))
    runsSel.value = new Set([...runsSel.value].filter((n) => names.has(n)))
  } catch (e) {
    flash(t('web.ml_results.load_fail', { msg: e instanceof Error ? e.message : String(e) }))
  }
}

const runsFiltered = computed(() => {
  const q = runsQuery.value.trim().toLowerCase()
  if (!q) return runs.value
  return runs.value.filter((r) =>
    r.name.toLowerCase().includes(q) ||
    String(r.gse ?? '').toLowerCase().includes(q) ||
    String(r.source ?? '').toLowerCase().includes(q))
})
const runsPaged = computed(() =>
  runsFiltered.value.slice((runsPage.value - 1) * RUNS_PAGE, runsPage.value * RUNS_PAGE))
const runsNPages = computed(() => Math.max(1, Math.ceil(runsFiltered.value.length / RUNS_PAGE)))
function runsGoto(p: number): void {
  runsPage.value = Math.min(Math.max(1, p), runsNPages.value)
}
function runsToggle(name: string): void {
  const next = new Set(runsSel.value)
  if (next.has(name)) next.delete(name)
  else next.add(name)
  runsSel.value = next
}
const runsAllSelected = computed(() =>
  runs.value.length > 0 && runs.value.every((r) => runsSel.value.has(r.name)))
function runsToggleAll(): void {
  runsSel.value = runsAllSelected.value
    ? new Set() : new Set(runs.value.map((r) => r.name))
}
function runsPaths(names: Iterable<string>): string[] {
  return [...names].map((n) => runs.value.find((r) => r.name === n)?.path ?? '').filter(Boolean)
}
async function runsZip(): Promise<void> {
  const paths = runsPaths(runsSel.value)
  if (!paths.length) return
  zipBusy.value = true
  try {
    const blob = await api.mlRunsZip(paths)
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `ml_results_${new Date().toISOString().slice(0, 10)}.zip`
    a.click()
    URL.revokeObjectURL(url)
  } catch {
    flash(t('web.common.zip_fail'))
  } finally {
    zipBusy.value = false
  }
}
async function runsDelete(): Promise<void> {
  if (!runsConfirm.value) {
    runsConfirm.value = true
    window.clearTimeout(runsConfirmTimer)
    runsConfirmTimer = window.setTimeout(() => (runsConfirm.value = false), 4000)
    return
  }
  runsConfirm.value = false
  const paths = runsPaths(runsSel.value)
  if (!paths.length) return
  runsDeleteBusy.value = true
  runsDeleteErr.value = ''
  try {
    const r = await api.mlRunsDelete({ paths })
    flash(t('web.ml_results.delete_done', { n: r.deleted.length }))
    runsSel.value = new Set()
    await loadRuns()
  } catch (e) {
    runsDeleteErr.value = e instanceof Error ? e.message : String(e)
  } finally {
    runsDeleteBusy.value = false
  }
}
async function runsCleanup(): Promise<void> {
  cleanupBusy.value = true
  cleanupErr.value = ''
  try {
    const r = await api.mlRunsCleanup({ keep: cleanupKeep.value })
    flash(t('web.ml_results.keepn_done', { n: r.deleted }))
    await loadRuns()
  } catch (e) {
    cleanupErr.value = e instanceof Error ? e.message : String(e)
  } finally {
    cleanupBusy.value = false
  }
}
const runsTotalMb = computed(() =>
  (runs.value.reduce((s, r) => s + Number(r.size_bytes ?? 0), 0) / 1048576).toFixed(0))
const fmtMb = (n: number): string => `${(Number(n) / 1048576).toFixed(1)} MB`
const fmtTime = (mtime: number): string =>
  new Date(Number(mtime) * 1000).toISOString().slice(0, 16).replace('T', ' ')

onMounted(() => void loadRuns())
</script>

<template>
  <p class="kicker s1">{{ t('web.ml_results.kicker') }}</p>
  <h1 class="title s2">{{ t('web.ml_results.title') }}</h1>
  <p class="sub s3">{{ t('web.ml_results.sub') }}</p>

  <div class="panel wide s4">
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <p v-if="runsLoaded && !runs.length" class="empty">{{ t('web.ml_results.empty_runs') }}</p>

    <div class="controls">
      <label class="search">
        <input v-model="runsQuery" spellcheck="false" :placeholder="t('web.ml_results.search_ph')"
               @input="runsGoto(1)" />
        <i class="line"></i>
      </label>
      <span class="count">{{ t('web.ml_results.runs_meta', { n: runsFiltered.length, mb: runsTotalMb, k: runsSel.size }) }}</span>
      <button class="go" :title="t('web.common.refresh')" @click="loadRuns">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <div class="bar" :class="{ dim: !runsSel.size }">
      <button class="bar-btn" @click="runsToggleAll">{{ runsAllSelected ? t('web.common.unselect_all') : t('web.common.all') }}</button>
      <span class="bar-label">{{ runsSel.size ? t('web.common.selected', { n: runsSel.size }) : t('web.common.batch_hint') }}</span>
      <button class="bar-btn" :disabled="!runsSel.size || zipBusy" @click="runsZip">
        {{ zipBusy ? t('web.ml_results.zipping') : t('web.common.download_selected') }}
      </button>
      <button class="bar-btn danger" :class="{ arm: runsConfirm }" :disabled="!runsSel.size || runsDeleteBusy"
              @click="runsDelete">
        {{ runsConfirm ? t('web.common.confirm_delete', { n: runsSel.size }) : t('web.common.delete_selected') }}
      </button>
      <span class="bar-label">{{ t('web.ml_results.keep_label') }}</span>
      <Stepper v-model="cleanupKeep" :min="1" :max="20" />
      <button class="bar-btn" :disabled="cleanupBusy" @click="runsCleanup">
        {{ cleanupBusy ? t('web.ml_results.cleaning') : t('web.ml_results.keepn_btn') }}
      </button>
    </div>
    <p v-if="runsDeleteErr" class="stale">{{ t('web.ml_results.delete_fail', { msg: runsDeleteErr }) }}</p>
    <p v-if="cleanupErr" class="stale">{{ t('web.ml_results.cleanup_fail', { msg: cleanupErr }) }}</p>

    <ul class="files">
      <li v-if="runsLoaded && !runsFiltered.length" class="empty">{{ t('web.common.no_match') }}</li>
      <li v-for="r in runsPaged" :key="r.name" class="file">
        <div class="row">
          <label class="sel">
            <input type="checkbox" :checked="runsSel.has(r.name)" @change="runsToggle(r.name)" />
          </label>
          <span class="fname" :title="r.name">{{ r.name
          }}<i v-if="r.stale?.length" class="stale-tag">{{ t('web.ml_results.stale_tag') }}</i></span>
          <span class="fmeta">
            {{ r.gse || r.source || '-' }} · {{ t('web.common.samples_n', { n: r.n_samples ?? '?' }) }} ·
            AUC {{ r.auc ?? 'N/A' }} · shortlist {{ r.shortlist_count ?? '?' }} ·
            {{ fmtMb(r.size_bytes) }} · {{ fmtTime(r.mtime) }}
          </span>
        </div>
      </li>
    </ul>

    <div v-if="runsNPages > 1" class="pager">
      <button class="pager-btn" :disabled="runsPage <= 1" @click="runsGoto(runsPage - 1)">‹</button>
      <span class="pager-num">{{ runsPage }} / {{ runsNPages }}</span>
      <button class="pager-btn" :disabled="runsPage >= runsNPages" @click="runsGoto(runsPage + 1)">›</button>
    </div>
  </div>
</template>

<style scoped>
.stale-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.2em;
  color: #b3261e;
  border: 1px solid #b3261e;
  padding: 0 5px;
  margin-left: 10px;
  vertical-align: 1.5px;
}
.stale {
  font-size: 12px;
  letter-spacing: 0.05em;
  color: #b3261e;
}
</style>
