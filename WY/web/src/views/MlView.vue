<script setup lang="ts">
/** GEO 表达谱 + ML 共识筛选(阶段 3c,Expression_ML_Screening 的 TS 迁移)。
 * 忠实迁移原页面四区布局:① 输入(GEO 检索/直填/本地矩阵) ② 配置(分组正则/排序表/参数)
 * ③ 运行(七阶段 job) ④ 结果(指标/ROC/共识/漏斗/shortlist + 近 3 次摘要)。
 * 末段并入 ML_Results_Information(结果管理):run 清单/搜索/删除/ZIP/keep-N 清理。
 * 全链编排唯一定义在 core.ml_pipeline.run_screen_chain(CLI/页面/API 共用)。 */
import { computed, onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl, type MlHistoryRun, type MlOverviewResult,
         type MlProbeResult, type MlRunResult, type Row } from '../api'
import Stepper from '../components/Stepper.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import { fmtSize } from '../utils/format'

const { t } = useI18n()

/* ---------------- ① 输入 ---------------- */

const tab = ref<'search' | 'direct' | 'local'>('search')
const diseases = ref<string[]>([])
const disease = ref('')
const limit = ref(15)
const searchResults = ref<Row[]>([])
const searchDisease = ref('')
const searchBusy = ref(false)
const gse = ref('')
const probe = ref<MlProbeResult | null>(null)
const probeBusy = ref(false)
let probedGse = ''
const exprPath = ref('')
const labelsPath = ref('')
const localCase = ref('')

async function loadDiseases(): Promise<void> {
  try {
    diseases.value = (await api.mlDiseases()).diseases
  } catch { /* 非阻断 */ }
}

async function search(): Promise<void> {
  if (!disease.value.trim()) { flash(t('web.ml.need_disease')); return }
  searchBusy.value = true
  try {
    const r = await api.mlSearch(disease.value.trim(), limit.value)
    searchResults.value = r.results as Row[]
    searchDisease.value = r.disease
    gse.value = ''
    probe.value = null
    probedGse = ''
  } catch (e) {
    flash(t('web.ml.search_fail', { msg: e instanceof Error ? e.message : String(e) }))
  } finally {
    searchBusy.value = false
  }
}

const gseOptions = computed<string[]>(() =>
  searchResults.value.map((r) => String(r.gse ?? '').trim()).filter(Boolean))
const searchStale = computed(() =>
  searchResults.value.length > 0 &&
  disease.value.trim().toLowerCase() !== searchDisease.value.toLowerCase())

function searchPlatformFor(g: string): { platform: string; type: string } {
  const row = searchResults.value.find((r) => String(r.gse ?? '').toUpperCase() === g)
  return { platform: String(row?.platform ?? ''), type: String(row?.type ?? '') }
}

/* GSE 确定后预检(每 GSE 一次) */
async function probeGse(): Promise<void> {
  const g = gse.value.trim().toUpperCase()
  if (!g || probedGse === g) return
  probedGse = g
  probeBusy.value = true
  try {
    probe.value = await api.mlProbe(g, searchPlatformFor(g))
  } catch (e) {
    probe.value = { err: e instanceof Error ? e.message : String(e) }
  } finally {
    probeBusy.value = false
  }
}

function onGseChange(): void {
  probedGse = ''
  probe.value = null
  void probeGse()
}

const useLocal = computed(() => !!exprPath.value.trim() && !!labelsPath.value.trim())
const localPartial = computed(() => !!exprPath.value.trim() !== !!labelsPath.value.trim())
const hasGse = computed(() => !!gse.value.trim())
const sourceConflict = computed(() => useLocal.value && hasGse.value)

/* 预检展示与可用性锁定 */
const probeR = computed(() => probe.value?.r ?? null)
const probeErr = computed(() => probe.value?.err ?? '')
const probeAssess = computed(() => probe.value?.assess ?? null)
const unusable = computed(() => {
  const p = probe.value
  const a = probeAssess.value
  return !!(hasGse.value && a && p?.gse === gse.value.trim().toUpperCase() &&
    ((a.matrix_empty && !a.counts_available) ||
     (!a.symbol_level && a.annot === 'missing')))
})

const fmtEta = (s: number): string =>
  s < 90 ? `${Math.max(1, Math.round(s))} ${t('web.ml.unit_sec')}`
    : s < 3600 ? `${Math.max(1, Math.round(s / 60))} ${t('web.ml.unit_min')}`
      : `${Math.floor(s / 3600)} ${t('web.ml.unit_hour')} ${Math.round(s % 3600 / 60)} ${t('web.ml.unit_min_short')}`

/* ---------------- ② 配置 ---------------- */

const casePat = ref('')
const ctrlPat = ref('')
const overviewJob = useJob((result) => {
  overview.value = result as unknown as MlOverviewResult
  flash(t('web.ml.overview_done_flash'))
})
const overview = ref<MlOverviewResult | null>(null)

async function loadOverview(): Promise<void> {
  const g = gse.value.trim().toUpperCase()
  if (!g) { flash(t('web.ml.need_gse')); return }
  await overviewJob.run(api.mlOverview(g), { doneText: t('web.ml.overview_done_text') })
}

/* 排序 CSV 候选 */
const rankedFiles = ref<{ label: string; base_disease: string; drug: string }[]>([])
const rankedSel = ref<Set<string>>(new Set())
const drugFilter = ref('')

const rankedShown = computed(() => {
  const ctx = disease.value.trim().toLowerCase().replace(/\s*rwr\s*\d+\s*$/i, '')
  const fd = drugFilter.value.trim().toLowerCase()
  return rankedFiles.value.filter((rf) =>
    (!ctx || ctx === rf.base_disease.toLowerCase() || ctx.includes(rf.base_disease.toLowerCase()))
    && (!fd || rf.drug.toLowerCase().includes(fd)))
})

async function loadRanked(): Promise<void> {
  try {
    const r = await api.mlRanked({ disease: disease.value, drug: drugFilter.value })
    rankedFiles.value = r.files
    rankedSel.value = new Set([...rankedSel.value].filter(
      (l) => r.files.some((f) => f.label === l)))
  } catch { /* 非阻断 */ }
}

function toggleRanked(label: string): void {
  const next = new Set(rankedSel.value)
  if (next.has(label)) next.delete(label)
  else next.add(label)
  rankedSel.value = next
}

/** 全选/清空只作用于当前过滤后可见的条目(滚动列表配套操作)。 */
function rankedSelectAll(on: boolean): void {
  const next = new Set(rankedSel.value)
  for (const f of rankedShown.value) {
    if (on) next.add(f.label)
    else next.delete(f.label)
  }
  rankedSel.value = next
}

/* 参数 */
const topK = ref(50)
const prefilterN = ref(500)
const bootstrapN = ref(30)

/* ---------------- ③ 运行 ---------------- */

const resultData = ref<MlRunResult | null>(null)
const fallback = ref<{ gse: string; date: string } | null>(null)
const runJob = useJob((result) => {
  resultData.value = result as unknown as MlRunResult
  const r = resultData.value
  if (r.geo_fallback) {
    fallback.value = { gse: String(r.gse), date: String(r.geo_cache_date ?? '') }
  } else {
    fallback.value = null
  }
  flash(t('web.ml.run_done_flash', { gse: r.gse, s: r.elapsed }))
  void loadHistory()
})

/* 切换数据上下文时清掉上一轮的展示,避免④区残留旧结果:
 * - 换 GSE:结果/总览/回退提示清空,分组正则一并清(样本命名规则变了),
 *   但三个调节器参数保留——那是用户的分析偏好,与数据集无关;
 * - 换疾病名 / 换输入方式(tab):只清结果区。 */
watch(gse, () => {
  resultData.value = null
  fallback.value = null
  overview.value = null
  casePat.value = ''
  ctrlPat.value = ''
})
watch(disease, () => {
  resultData.value = null
  fallback.value = null
})
watch(tab, () => {
  resultData.value = null
  fallback.value = null
})

/* 整包下载:当前 run 目录打成一个 zip(单 run 走 make_dir_zip) */
const zipBusy = ref(false)
const runDirName = computed(() => (resultData.value?.out_dir ?? '').split(/[\\/]/).pop() ?? '')
async function zipRun(): Promise<void> {
  if (!resultData.value?.out_dir) return
  zipBusy.value = true
  try {
    const blob = await api.mlRunsZip([resultData.value.out_dir])
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${runDirName.value}.zip`
    a.click()
    URL.revokeObjectURL(url)
  } catch {
    flash(t('web.common.zip_fail'))
  } finally {
    zipBusy.value = false
  }
}

const canRun = computed(() =>
  ((hasGse.value && casePat.value.trim() && ctrlPat.value.trim()) || useLocal.value)
  && !sourceConflict.value && !unusable.value)

async function run(): Promise<void> {
  if (sourceConflict.value) { flash(t('web.ml.conflict_flash')); return }
  if (localPartial.value) { flash(t('web.ml.partial_flash')); return }
  await runJob.run(api.mlRun({
    mode: useLocal.value ? 'local' : 'gse',
    gse: gse.value.trim().toUpperCase(),
    case_pattern: casePat.value, control_pattern: ctrlPat.value,
    expr: exprPath.value, labels: labelsPath.value, case_label: localCase.value,
    ranked_labels: [...rankedSel.value],
    top_k: topK.value, prefilter_top: prefilterN.value, bootstrap_rounds: bootstrapN.value
  }), { doneText: t('web.ml.run_done_text') })
}

const forceJob = useJob(() => {
  fallback.value = null
  probedGse = ''
  probe.value = null
  flash(t('web.ml.force_done_flash'))
})

async function forceRedownload(): Promise<void> {
  const g = fallback.value?.gse ?? gse.value.trim().toUpperCase()
  if (!g) return
  await forceJob.run(api.mlForceRefresh(g), { doneText: t('web.ml.force_done_text') })
}

/* ---------------- ④ 结果 ---------------- */

const meta = computed(() => resultData.value?.meta ?? null)
const consensusRows = computed<Row[]>(() => resultData.value?.consensus ?? [])
const consensusTotal = computed(() => resultData.value?.consensus_total ?? 0)
const funnelRows = computed<Row[]>(() => resultData.value?.funnel ?? [])
const shortlistRows = computed<Row[]>(() => resultData.value?.shortlist ?? [])
const shortlistCount = computed(() => resultData.value?.shortlist_count ?? '—')
const artifactsMap = computed<Record<string, string>>(() => resultData.value?.artifacts ?? {})

/* ---------------- 下载区(run 目录全部落盘文件,按说明归类) ---------------- */

interface DlItem { file: string; rel: string; desc: string; main: boolean }

/** 产物归类:页面展示对应的四项为主成果;DEG 全表/共识曲线/元数据等为补充产物(默认折叠)。
 * 未识别的新产物文件名即说明,进补充区——后端新增产物不需要改前端就能出现。 */
function classifyArtifact(file: string): { desc: string; main: boolean } | null {
  if (/^ml_roc_.*\.png$/.test(file)) return { desc: t('web.ml.dl_roc'), main: true }
  if (file === 'ml_consensus_scores.csv') return { desc: t('web.ml.dl_consensus'), main: true }
  if (file === 'shortlist.csv') return { desc: t('web.ml.dl_shortlist'), main: true }
  if (file === 'funnel.csv') return { desc: t('web.ml.dl_funnel'), main: true }
  if (file === 'deg_results.csv') return { desc: t('web.ml.dl_deg'), main: false }
  if (file === 'consensus_curve.png') return { desc: t('web.ml.dl_curve'), main: false }
  if (file === 'summary.json') return { desc: t('web.ml.dl_summary'), main: false }
  return null
}

const mainDownloads = computed<DlItem[]>(() => {
  const out: DlItem[] = []
  for (const [file, rel] of Object.entries(artifactsMap.value)) {
    const c = classifyArtifact(file)
    if (c?.main) out.push({ file, rel, desc: c.desc, main: true })
  }
  return out
})
const extraDownloads = computed<DlItem[]>(() => {
  const out: DlItem[] = []
  for (const [file, rel] of Object.entries(artifactsMap.value)) {
    const c = classifyArtifact(file)
    if (!c?.main) out.push({ file, rel, desc: c?.desc ?? file, main: false })
  }
  return out
})

/* ---------------- 历史(近 3 次) ---------------- */

const history = ref<MlHistoryRun[]>([])
const histPick = ref('')
const histRow = computed(() => history.value.find((r) => r.name === histPick.value) ?? null)

async function loadHistory(): Promise<void> {
  try {
    history.value = (await api.mlHistory()).runs
    if (!histPick.value && history.value.length) histPick.value = history.value[0].name
  } catch { /* 非阻断 */ }
}

/* ---------------- 通用 ---------------- */

const { notice, flash } = useNotice(4200)

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

onMounted(() => {
  void loadDiseases()
  void loadRanked()
  void loadHistory()
})
</script>

<template>
  <p class="kicker s1">{{ t('web.ml.kicker') }}</p>
  <h1 class="title s2">{{ t('web.ml.title') }}</h1>
  <p class="sub s3">{{ t('web.ml.sub') }}</p>

  <div class="panel wide s4">
    <div class="xlinks">
      <a class="bar-btn" href="#/ml-results">{{ t('web.nav.results_link') }}</a>
    </div>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <!-- ============ ① 输入区 ============ -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.ml.sec_input') }}</h3>
      <div class="lv-group">
        <button class="lv" :class="{ on: tab === 'search' }" @click="tab = 'search'">{{ t('web.ml.tab_search') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: tab === 'direct' }" @click="tab = 'direct'">{{ t('web.ml.tab_direct') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: tab === 'local' }" @click="tab = 'local'">{{ t('web.ml.tab_local') }}</button>
      </div>

      <template v-if="tab === 'search'">
        <div class="set-row">
          <label class="set-label">{{ t('web.ml.disease_label') }}</label>
          <div class="field wide">
            <input v-model="disease" list="ml-diseases" spellcheck="false"
                   :placeholder="t('web.ml.disease_ph')" @change="loadRanked" />
            <i class="line"></i>
            <datalist id="ml-diseases">
              <option v-for="d in diseases" :key="d" :value="d" />
            </datalist>
          </div>
          <span class="ctl-label">{{ t('web.ml.limit_label') }}</span>
          <Stepper v-model="limit" :min="5" :max="50" />
          <button class="bar-btn strong" :disabled="searchBusy" @click="search">
            {{ searchBusy ? t('web.ml.searching') : t('web.ml.search_btn') }}
          </button>
        </div>
        <p v-if="diseases.length" class="hint-line">{{ t('web.ml.diseases_hint', { n: diseases.length }) }}</p>
        <div v-if="searchResults.length" class="tbl-wrap">
          <p v-if="searchStale" class="stale">{{ t('web.ml.search_stale', { d: searchDisease }) }}</p>
          <table class="tbl">
            <thead><tr><th v-for="k in Object.keys(searchResults[0])" :key="k">{{ k }}</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in searchResults" :key="i">
                <td v-for="(v, k) in r" :key="k">{{ v }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p v-else-if="searchDisease" class="hint-line">{{ t('web.ml.no_results') }}</p>
      </template>
      <p v-else-if="tab === 'direct'" class="hint-line">{{ t('web.ml.direct_hint') }}</p>

      <template v-if="tab === 'local'">
        <p class="hint-line">{{ t('web.ml.local_hint') }}</p>
        <div class="set-row">
          <label class="set-label">{{ t('web.ml.expr_label') }}</label>
          <div class="field wide"><input v-model="exprPath" spellcheck="false"
                                         placeholder="D:\data\GSE12345_expr.csv" /><i class="line"></i></div>
        </div>
        <div class="set-row">
          <label class="set-label">{{ t('web.ml.labels_label') }}</label>
          <div class="field wide"><input v-model="labelsPath" spellcheck="false"
                                         placeholder="D:\data\labels.csv" /><i class="line"></i></div>
        </div>
        <div class="set-row">
          <label class="set-label">{{ t('web.ml.case_label_value') }}</label>
          <div class="field"><input v-model="localCase" spellcheck="false"
                                    :placeholder="t('web.ml.case_ph')" /><i class="line"></i></div>
        </div>
      </template>

      <!-- GSE 确认(检索/直填共用) -->
      <div class="set-row" v-if="tab !== 'local'">
        <label class="set-label">{{ t('web.ml.gse_label') }}</label>
        <select v-if="gseOptions.length" class="sel" v-model="gse" @change="onGseChange">
          <option value="" disabled>{{ t('web.ml.sel_gse_ph') }}</option>
          <option v-for="g in gseOptions" :key="g" :value="g">{{ g }}</option>
        </select>
        <div v-else class="field">
          <input v-model="gse" placeholder="GSE63067" spellcheck="false"
                 @change="onGseChange" /><i class="line"></i>
        </div>
        <span v-if="gse" class="state on">{{ t('web.ml.confirmed', { gse }) }}</span>
      </div>

      <!-- 预检展示 -->
      <template v-if="tab !== 'local' && gse && (probeBusy || probe)">
        <p v-if="probeBusy" class="hint-line">{{ t('web.ml.probing') }}</p>
        <template v-else>
          <p v-if="probeR" class="hint-line">
            {{ t('web.ml.probe_meta', { size: fmtSize(probeR.size_bytes),
                                        speed: (probeR.speed_bps / 1024).toFixed(0),
                                        eta: fmtEta(Math.min(probeR.size_bytes / Math.max(probeR.speed_bps, 1), 5400)) }) }}
          </p>
          <p v-if="probeR && (!probeR.reliable || probeR.size_bytes / Math.max(probeR.speed_bps, 1) > 300)"
             class="stale">{{ t('web.ml.probe_slow') }}</p>
          <p v-if="probeErr" class="hint-line">{{ t('web.ml.probe_fail', { msg: probeErr }) }}</p>
          <template v-if="probeAssess">
            <p v-if="probeAssess.matrix_empty && probeAssess.counts_available" class="hint-line good">
              {{ t('web.ml.counts_fallback') }}
            </p>
            <p v-else-if="probeAssess.matrix_empty" class="stale">{{ t('web.ml.empty_unusable') }}</p>
            <p v-else-if="!probeAssess.symbol_level && probeAssess.annot === 'missing'"
               class="stale">{{ t('web.ml.probe_annot_missing') }}</p>
          </template>
          <p v-else-if="probe?.annot_missing" class="stale">
            {{ t('web.ml.annot_unreachable', { gpl: probe?.gpl }) }}
          </p>
        </template>
      </template>
    </section>

    <!-- ============ ② 配置区 ============ -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.ml.sec_config') }}</h3>
      <p v-if="unusable" class="stale">{{ t('web.ml.config_locked') }}</p>
      <div class="set-row">
        <label class="set-label">{{ t('web.ml.case_pat_label') }}</label>
        <div class="field wide"><input v-model="casePat" :disabled="unusable"
                                       placeholder="steatosis|non-alcoholic" spellcheck="false" /><i class="line"></i></div>
      </div>
      <div class="set-row">
        <label class="set-label">{{ t('web.ml.ctrl_pat_label') }}</label>
        <div class="field wide"><input v-model="ctrlPat" :disabled="unusable"
                                       placeholder="healthy" spellcheck="false" /><i class="line"></i></div>
        <button class="bar-btn" :disabled="!gse || overviewJob.running.value || unusable"
                @click="loadOverview">
          {{ overviewJob.running.value ? t('web.ml.overview_running') : t('web.ml.overview_btn') }}
        </button>
      </div>
      <div v-if="overviewJob.running.value" class="progress">
        <div class="progress-track"><div class="progress-fill indeterminate"></div></div>
        <p class="progress-text">{{ overviewJob.text.value || t('web.ml.overview_dl') }}</p>
      </div>
      <p v-if="overviewJob.error.value" class="stale">{{ t('web.ml.overview_fail', { msg: overviewJob.error.value }) }}</p>

      <template v-if="overview">
        <div class="tbl-wrap">
          <table class="tbl">
            <thead><tr><th v-for="k in Object.keys(overview.overview[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in overview.overview" :key="i">
                <td v-for="(v, k) in r" :key="k">{{ v }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p v-if="overview.issues?.length && !overview.counts_rows?.length" class="stale">
          {{ t('web.ml.issues', { issues: overview.issues.join('; ') }) }}
        </p>
        <template v-if="overview.counts_rows?.length">
          <p class="hint-line good">{{ t('web.ml.counts_hint', { n: overview.counts_rows.length }) }}</p>
          <div class="tbl-wrap" style="max-height:260px">
            <table class="tbl">
              <thead><tr><th v-for="k in Object.keys(overview.counts_rows[0])" :key="k">{{ k }}</th></tr></thead>
              <tbody>
                <tr v-for="(r, i) in overview.counts_rows" :key="i">
                  <td v-for="(v, k) in r" :key="k">{{ v }}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </template>
        <template v-if="overview.assess && !overview.assess.matrix_empty">
          <p class="hint-line good">
            {{ t('web.ml.matrix_meta', {
              r: overview.assess.shape?.[0] ?? '—', c: overview.assess.shape?.[1] ?? '—',
              level: overview.assess.symbol_level ? t('web.ml.id_symbol') : t('web.ml.id_probe'),
              annot: overview.assess.symbol_level ? ''
                : (overview.assess.annot === 'cached' || overview.assess.annot === 'ok'
                  ? t('web.ml.annot_ok', { gpl: overview.assess.gpl })
                  : overview.assess.annot === 'missing' ? t('web.ml.annot_missing', { gpl: overview.assess.gpl })
                    : t('web.ml.annot_unknown', { gpl: overview.assess.gpl }))
            }) }}
          </p>
          <p v-if="overview.assess.group_fields?.length" class="hint-line">
            {{ t('web.ml.group_fields', { fields: overview.assess.group_fields.slice(0, 3).join('; ') }) }}
          </p>
          <p v-if="overview.assess.n_samples != null && overview.assess.n_samples < 12" class="stale">
            {{ t('web.ml.few_samples', { n: overview.assess.n_samples }) }}
          </p>
        </template>
      </template>

      <!-- 排序 CSV -->
      <div class="set-row" style="margin-top:16px">
        <label class="set-label">{{ t('web.ml.ranked_label') }}</label>
        <div class="field"><input v-model="drugFilter" spellcheck="false"
                                  :placeholder="t('web.ml.drug_filter_ph')" @change="loadRanked" /><i class="line"></i></div>
        <span class="ctl-meta">{{ t('web.ml.ranked_meta', { n: rankedShown.length, m: rankedSel.size }) }}</span>
        <template v-if="rankedShown.length">
          <button type="button" class="bar-btn" :disabled="unusable"
                  @click="rankedSelectAll(true)">{{ t('web.common.all') }}</button>
          <button type="button" class="bar-btn" :disabled="unusable"
                  @click="rankedSelectAll(false)">{{ t('web.ml.ranked_clear') }}</button>
        </template>
      </div>
      <!-- 定高滚动列表:候选多时不再折行挤压,一行一项 -->
      <div v-if="rankedShown.length" class="ranked-list">
        <label v-for="f in rankedShown" :key="f.label" class="ranked-item">
          <input type="checkbox" :checked="rankedSel.has(f.label)" @change="toggleRanked(f.label)"
                 :disabled="unusable" />
          <span class="ranked-name" :title="f.label">{{ f.label }}</span>
        </label>
      </div>
      <p v-else-if="disease" class="hint-line">{{ t('web.ml.ranked_none_disease') }}</p>
      <p v-else class="hint-line">{{ t('web.ml.ranked_hint') }}</p>

      <div class="opt-grid">
        <label class="opt">
          <span class="opt-label">{{ t('web.ml.opt_topk') }}</span>
          <Stepper v-model="topK" :min="5" :max="100" />
        </label>
        <label class="opt">
          <span class="opt-label">{{ t('web.ml.opt_prefilter') }}</span>
          <Stepper v-model="prefilterN" :min="50" :max="5000" :step="50" />
        </label>
        <label class="opt">
          <span class="opt-label">{{ t('web.ml.opt_bootstrap') }}</span>
          <Stepper v-model="bootstrapN" :min="5" :max="200" :step="5" />
        </label>
      </div>
      <p class="hint-line">{{ t('web.ml.perf_hint') }}</p>
    </section>

    <!-- ============ ③ 运行区 ============ -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.ml.sec_run') }}</h3>
      <p v-if="unusable" class="stale">{{ t('web.ml.run_locked') }}</p>
      <p v-else-if="sourceConflict" class="stale">{{ t('web.ml.conflict') }}</p>
      <p v-else-if="localPartial" class="stale">{{ t('web.ml.local_partial') }}</p>
      <button class="build" :disabled="!canRun || runJob.running.value || overviewJob.running.value"
              @click="run">
        <svg viewBox="0 0 24 14" width="26" height="15">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <span class="build-hint">{{ t('web.ml.run_hint') }}</span>

      <div v-if="runJob.running.value" class="prog">
        <div class="progress">
          <div class="progress-track">
            <div class="progress-fill"
                 :style="{ width: `${((runJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
          </div>
        </div>
        <ul class="stages">
          <li v-for="k in STAGE_ORDER" :key="k" class="stage" :class="stageStates[k]">
            <i class="stage-icon">{{ STATE_ICON[stageStates[k]] }}</i>{{ stageLabels[k] }}
          </li>
        </ul>
        <p v-if="runJob.text.value" class="hint-line">{{ runJob.text.value }}</p>
      </div>
      <p v-if="runJob.error.value" class="stale">{{ t('web.ml.run_fail', { msg: runJob.error.value }) }}</p>

      <div v-if="fallback" class="set-row" style="margin-top:14px">
        <span class="stale">{{ t('web.ml.fallback_warn', { date: fallback.date }) }}</span>
        <button class="bar-btn" :disabled="forceJob.running.value" @click="forceRedownload">
          {{ forceJob.running.value ? t('web.ml.forcing') : t('web.ml.force_btn') }}
        </button>
      </div>
      <p v-if="forceJob.error.value" class="stale">{{ t('web.ml.force_fail', { msg: forceJob.error.value }) }}</p>
    </section>

    <!-- ============ ④ 结果区 ============ -->
    <section class="sec" v-if="resultData">
      <h3 class="sec-title">{{ t('web.ml.sec_result', { gse: resultData.gse, s: resultData.elapsed }) }}</h3>
      <div class="metrics">
        <div class="metric"><b>{{ meta?.n_samples }}</b><span>{{ t('web.ml.m_samples') }}</span></div>
        <div class="metric"><b>{{ meta?.n_genes }}</b><span>{{ t('web.ml.m_genes') }}</span></div>
        <div class="metric"><b>{{ meta?.cv_auc ?? 'N/A' }}</b><span>CV-AUC</span></div>
        <div class="metric"><b>{{ meta?.algo_spearman_mean }}</b><span>{{ t('web.ml.m_spearman') }}</span></div>
      </div>

      <template v-if="resultData.roc_png_rel">
        <p class="hint-line" style="margin-top:16px">{{ t('web.ml.roc_hint') }}</p>
        <img class="fig-img" :src="artifactsUrl(resultData.roc_png_rel)" alt="ROC" />
      </template>

      <h4 class="sub-sec">{{ t('web.ml.consensus_head', { n: Math.min(20, consensusRows.length), total: consensusTotal }) }}</h4>
      <div class="tbl-wrap">
        <table class="tbl">
          <thead><tr><th v-for="k in Object.keys(consensusRows[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
          <tbody>
            <tr v-for="(r, i) in consensusRows.slice(0, 20)" :key="i">
              <td v-for="(v, k) in r" :key="k">{{ v }}</td>
            </tr>
          </tbody>
        </table>
      </div>

      <h4 class="sub-sec">{{ t('web.ml.funnel_head') }}</h4>
      <div class="tbl-wrap" style="max-height:240px">
        <table class="tbl">
          <thead><tr><th v-for="k in Object.keys(funnelRows[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
          <tbody><tr v-for="(r, i) in funnelRows" :key="i"><td v-for="(v, k) in r" :key="k">{{ v }}</td></tr></tbody>
        </table>
      </div>

      <h4 class="sub-sec">{{ t('web.ml.shortlist_head', { n: shortlistCount }) }}</h4>
      <div class="tbl-wrap" style="max-height:320px">
        <table class="tbl">
          <thead><tr><th v-for="k in Object.keys(shortlistRows[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
          <tbody><tr v-for="(r, i) in shortlistRows" :key="i"><td v-for="(v, k) in r" :key="k">{{ v }}</td></tr></tbody>
        </table>
      </div>

      <!-- 下载区:整包一键带走;主成果带说明(与分析页同款排版),补充产物默认折叠 -->
      <h4 class="sub-sec">{{ t('web.ml.dl_title') }}</h4>
      <button class="bar-btn" :disabled="zipBusy" @click="zipRun">
        {{ zipBusy ? t('web.ml_results.zipping') : t('web.ml.dl_zip', { name: `${runDirName}.zip` }) }}
      </button>
      <ul v-if="mainDownloads.length" class="dl-list">
        <li v-for="d in mainDownloads" :key="d.rel" class="dl-row">
          <div class="dl-main">
            <span class="dl-desc">{{ d.desc }}</span>
            <span class="dl-name">{{ d.file }}</span>
          </div>
          <a class="mini" :href="artifactsUrl(d.rel, true)">{{ t('web.common.download') }}</a>
        </li>
      </ul>
      <details v-if="extraDownloads.length" class="fold">
        <summary class="hint-line fold2">{{ t('web.ml.dl_more', { n: extraDownloads.length }) }}</summary>
        <ul class="dl-list">
          <li v-for="d in extraDownloads" :key="d.rel" class="dl-row">
            <div class="dl-main">
              <span class="dl-desc">{{ d.desc }}</span>
              <span class="dl-name">{{ d.file }}</span>
            </div>
            <a class="mini" :href="artifactsUrl(d.rel, true)">{{ t('web.common.download') }}</a>
          </li>
        </ul>
      </details>
    </section>
    <section class="sec" v-else>
      <h3 class="sec-title">{{ t('web.ml.sec_result_idle') }}</h3>
      <p class="hint-line">{{ t('web.ml.result_idle_hint') }}</p>
    </section>

    <!-- ============ 历史结果(近 3 次) ============ -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.ml.history_head') }}</h3>
      <p v-if="!history.length" class="hint-line">{{ t('web.ml.history_empty') }}</p>
      <div v-else class="set-row">
        <select class="sel inline" v-model="histPick">
          <option v-for="r in history" :key="r.name" :value="r.name">{{ r.name }}</option>
        </select>
        <span v-if="histRow" class="ctl-meta">
          {{ t('web.ml.hist_meta', { auc: histRow.auc ?? 'N/A', n: t('web.common.samples_n', { n: histRow.n_samples ?? '?' }),
                                     sl: histRow.shortlist_count ?? '?', ver: histRow.ml_version ?? '?' }) }}
        </span>
        <span v-for="s in histRow?.stale ?? []" :key="s.type" class="stale">
          {{ s.type === 'version' ? t('web.ml.stale_version', { old: s.old, new: s.new })
            : s.type === 'ranked_changed' ? t('web.ml.stale_ranked_changed', { name: s.name })
              : t('web.ml.stale_ranked_missing', { name: s.name }) }}
        </span>
      </div>
    </section>

    <p class="hint-line">
      {{ t('web.ml.manage_hint') }}
      <a class="bar-btn" href="#/ml-results">{{ t('web.ml.manage_link') }}</a>
    </p>
  </div>
</template>

<style scoped>
.state.on {
  color: #1a1a1a;
  font-weight: 500;
}
.build {
  width: 52px;
  height: 52px;
  border: none;
  border-radius: 50%;
  background: #1a1a1a;
  color: #fdfdfd;
  display: grid;
  place-items: center;
  cursor: pointer;
  transition: transform 0.3s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.build:hover:not(:disabled) {
  transform: scale(1.06);
}
.build:disabled {
  opacity: 0.4;
  cursor: default;
}
.build-hint {
  font-size: 11px;
  letter-spacing: 0.3em;
  color: #bcbcbc;
  margin-left: 14px;
}
.prog {
  margin-top: 18px;
}
.stages {
  list-style: none;
  margin: 10px 0 0;
  padding: 0;
  display: flex;
  flex-wrap: wrap;
  gap: 8px 20px;
}
.stage {
  font-size: 12px;
  letter-spacing: 0.06em;
  color: #c6c6c4;
}
.stage.running {
  color: #1a1a1a;
}
.stage.done {
  color: #8a8a88;
}
.stage-icon {
  margin-right: 6px;
  font-style: normal;
}
.fig-img {
  max-width: 420px;
  margin-top: 10px;
  border: 1px solid #f0f0ee;
  background: #ffffff;
}
.sub-sec {
  margin: 22px 0 8px;
  font-size: 13px;
  font-weight: 500;
  letter-spacing: 0.12em;
  color: #1a1a1a;
}
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
/* 平台排序表:定高滚动列表(发丝线边框,一行一项,候选多时滚动选) */
.ranked-list {
  margin-top: 10px;
  max-height: 176px;
  overflow-y: auto;
  border: 1px solid #f0f0ee;
}
.ranked-item {
  display: flex;
  align-items: baseline;
  gap: 10px;
  padding: 7px 12px;
  border-bottom: 1px solid #f4f4f2;
  font-size: 12px;
  letter-spacing: 0.04em;
  color: #4a4a48;
  cursor: pointer;
  user-select: none;
}
.ranked-item:last-child {
  border-bottom: none;
}
.ranked-item:hover {
  background: #f7f7f5;
}
.ranked-item input {
  accent-color: #1a1a1a;
  flex: none;
}
.ranked-name {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-variant-numeric: tabular-nums;
}
/* 排序表列表与下方参数调节器之间留出呼吸间距 */
.opt-grid {
  margin-top: 30px;
}
/* 补充产物折叠(下载区) */
.fold2 {
  margin-top: 12px;
  cursor: pointer;
  list-style: none;
  user-select: none;
}
.fold2::before {
  content: '› ';
  color: #c2c2c2;
}
.fold[open] .fold2::before {
  content: '⌄ ';
}
</style>
