<script setup lang="ts">
/** GEO 表达谱 + ML 共识筛选(阶段 3c,Expression_ML_Screening 的 TS 迁移)。
 * 忠实迁移原页面四区布局:① 输入(GEO 检索/直填/本地矩阵) ② 配置(分组正则/排序表/参数)
 * ③ 运行(七阶段 job) ④ 结果(指标/ROC/共识/漏斗/shortlist + 近 3 次摘要)。
 * 末段并入 ML_Results_Information(结果管理):run 清单/搜索/删除/ZIP/keep-N 清理。
 * 全链编排唯一定义在 core.ml_pipeline.run_screen_chain(CLI/页面/API 共用)。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { api, artifactsUrl, type Row } from '../api'
import Stepper from '../components/Stepper.vue'
import { useJob } from '../composables/useJob'

/* ---------------- ① 输入 ---------------- */

const tab = ref<'search' | 'direct' | 'local'>('search')
const diseases = ref<string[]>([])
const disease = ref('')
const limit = ref(15)
const searchResults = ref<Row[]>([])
const searchDisease = ref('')
const searchBusy = ref(false)
const gse = ref('')
const probe = ref<Record<string, any> | null>(null)
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
  if (!disease.value.trim()) { flash('请先填写疾病名称'); return }
  searchBusy.value = true
  try {
    const r = await api.mlSearch(disease.value.trim(), limit.value)
    searchResults.value = r.results as Row[]
    searchDisease.value = r.disease
    gse.value = ''
    probe.value = null
    probedGse = ''
  } catch (e) {
    flash('检索失败:' + (e instanceof Error ? e.message : String(e)))
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

const fmtSize = (n: number): string =>
  n < 1048576 ? `${(n / 1024).toFixed(0)} KB` : `${(n / 1048576).toFixed(1)} MB`
const fmtEta = (s: number): string =>
  s < 90 ? `${Math.max(1, Math.round(s))} 秒`
    : s < 3600 ? `${Math.max(1, Math.round(s / 60))} 分钟`
      : `${Math.floor(s / 3600)} 小时 ${Math.round(s % 3600 / 60)} 分`

/* ---------------- ② 配置 ---------------- */

const casePat = ref('')
const ctrlPat = ref('')
const overviewJob = useJob((result) => {
  overview.value = result as Record<string, any>
  flash('总览已生成(可用性判定已更新)')
})
const overview = ref<Record<string, any> | null>(null)

async function loadOverview(): Promise<void> {
  const g = gse.value.trim().toUpperCase()
  if (!g) { flash('先确认 GSE 编号'); return }
  await overviewJob.run(api.mlOverview(g), { doneText: '总览完成' })
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

/* 参数 */
const topK = ref(50)
const prefilterN = ref(500)
const bootstrapN = ref(30)

/* ---------------- ③ 运行 ---------------- */

const resultData = ref<Record<string, any> | null>(null)
const fallback = ref<{ gse: string; date: string } | null>(null)
const runJob = useJob((result) => {
  resultData.value = result as Record<string, any>
  const r = resultData.value
  if (r.geo_fallback) {
    fallback.value = { gse: String(r.gse), date: String(r.geo_cache_date ?? '') }
  } else {
    fallback.value = null
  }
  flash(`筛选完成:${r.gse} · 用时 ${r.elapsed}s`)
  void loadHistory()
})

const canRun = computed(() =>
  ((hasGse.value && casePat.value.trim() && ctrlPat.value.trim()) || useLocal.value)
  && !sourceConflict.value && !unusable.value)

async function run(): Promise<void> {
  if (sourceConflict.value) { flash('GSE 与本地矩阵二选一,不能同时提供'); return }
  if (localPartial.value) { flash('本地矩阵需要同时填写表达谱与标签路径'); return }
  await runJob.run(api.mlRun({
    mode: useLocal.value ? 'local' : 'gse',
    gse: gse.value.trim().toUpperCase(),
    case_pattern: casePat.value, control_pattern: ctrlPat.value,
    expr: exprPath.value, labels: labelsPath.value, case_label: localCase.value,
    ranked_labels: [...rankedSel.value],
    top_k: topK.value, prefilter_top: prefilterN.value, bootstrap_rounds: bootstrapN.value
  }), { doneText: '筛选完成' })
}

const forceJob = useJob(() => {
  fallback.value = null
  probedGse = ''
  probe.value = null
  flash('GEO 缓存已强制重下,预检将重新判定')
})

async function forceRedownload(): Promise<void> {
  const g = fallback.value?.gse ?? gse.value.trim().toUpperCase()
  if (!g) return
  await forceJob.run(api.mlForceRefresh(g), { doneText: '重下完成' })
}

/* ---------------- ④ 结果 ---------------- */

const meta = computed<Record<string, any> | null>(() => resultData.value?.meta ?? null)
const consensusRows = computed<Row[]>(() => resultData.value?.consensus ?? [])
const consensusTotal = computed(() => resultData.value?.consensus_total ?? 0)
const funnelRows = computed<Row[]>(() => resultData.value?.funnel ?? [])
const shortlistRows = computed<Row[]>(() => resultData.value?.shortlist ?? [])
const shortlistCount = computed(() => resultData.value?.shortlist_count ?? '—')
const artifactsMap = computed<Record<string, string>>(() => resultData.value?.artifacts ?? {})

/* ---------------- 历史(近 3 次) ---------------- */

const history = ref<Record<string, any>[]>([])
const histPick = ref('')
const histRow = computed(() => history.value.find((r) => r.name === histPick.value) ?? null)

async function loadHistory(): Promise<void> {
  try {
    history.value = (await api.mlHistory()).runs
    if (!histPick.value && history.value.length) histPick.value = history.value[0].name
  } catch { /* 非阻断 */ }
}

/* ---------------- 通用 ---------------- */

const notice = ref('')
let noticeTimer = 0
function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 4200)
}

const STAGE_ZH: Record<string, string> = {
  fetch: '获取数据集', group: '样本分组', preprocess: '预处理',
  deg: '差异表达', ml: 'ML 共识', shortlist: 'shortlist', save: '落盘'
}
const STATE_ICON: Record<string, string> = { pending: '○', running: '◐', done: '●' }
const stageOrder = Object.keys(STAGE_ZH)
const stageStates = computed<Record<string, string>>(() => {
  const map: Record<string, string> = {}
  for (const k of stageOrder) map[k] = 'pending'
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
  <p class="kicker s1">GEO&nbsp;EXPRESSION&nbsp;+&nbsp;ML&nbsp;CONSENSUS</p>
  <h1 class="title s2">表达谱 ML 共识筛选</h1>
  <p class="sub s3">检索&nbsp;—&nbsp;分组&nbsp;—&nbsp;DEG&nbsp;—&nbsp;ML&nbsp;—&nbsp;shortlist</p>

  <div class="panel wide s4">
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <!-- ============ ① 输入区 ============ -->
    <section class="sec">
      <h3 class="sec-title">① 数据来源</h3>
      <div class="lv-group">
        <button class="lv" :class="{ on: tab === 'search' }" @click="tab = 'search'">GEO 检索</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: tab === 'direct' }" @click="tab = 'direct'">直填 GSE</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: tab === 'local' }" @click="tab = 'local'">本地矩阵</button>
      </div>

      <template v-if="tab === 'search'">
        <div class="set-row">
          <label class="set-label">疾病名称</label>
          <div class="field wide">
            <input v-model="disease" list="ml-diseases" spellcheck="false"
                   placeholder="从已分析疾病中选择或直接输入" @change="loadRanked" />
            <i class="line"></i>
            <datalist id="ml-diseases">
              <option v-for="d in diseases" :key="d" :value="d" />
            </datalist>
          </div>
          <span class="ctl-label">条数</span>
          <Stepper v-model="limit" :min="5" :max="50" />
          <button class="bar-btn strong" :disabled="searchBusy" @click="search">
            {{ searchBusy ? '检索中…' : '检索 GEO' }}
          </button>
        </div>
        <p v-if="diseases.length" class="hint-line">已带出本机 {{ diseases.length }} 个分析过的疾病,可直接选择。</p>
        <div v-if="searchResults.length" class="tbl-wrap">
          <p v-if="searchStale" class="stale">当前疾病与检索时({{ searchDisease }})不一致,结果仅供参考。</p>
          <table class="tbl">
            <thead><tr><th v-for="k in Object.keys(searchResults[0])" :key="k">{{ k }}</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in searchResults" :key="i">
                <td v-for="(v, k) in r" :key="k">{{ v }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p v-else-if="searchDisease" class="hint-line">检索无结果。</p>
      </template>
      <p v-else-if="tab === 'direct'" class="hint-line">直接在下方输入 GSE 编号(如 GSE63067)。</p>

      <template v-if="tab === 'local'">
        <p class="hint-line">本地矩阵模式:表达谱 CSV(行=基因,列=样本)+ 标签 CSV(sample_id, group 两列)。两路都填时优先于 GSE。</p>
        <div class="set-row">
          <label class="set-label">表达谱 CSV</label>
          <div class="field wide"><input v-model="exprPath" spellcheck="false"
                                         placeholder="D:\data\GSE12345_expr.csv" /><i class="line"></i></div>
        </div>
        <div class="set-row">
          <label class="set-label">标签 CSV</label>
          <div class="field wide"><input v-model="labelsPath" spellcheck="false"
                                         placeholder="D:\data\labels.csv" /><i class="line"></i></div>
        </div>
        <div class="set-row">
          <label class="set-label">case 标签值</label>
          <div class="field"><input v-model="localCase" spellcheck="false"
                                    placeholder="如 disease(留空自动判定)" /><i class="line"></i></div>
        </div>
      </template>

      <!-- GSE 确认(检索/直填共用) -->
      <div class="set-row" v-if="tab !== 'local'">
        <label class="set-label">GSE 编号</label>
        <select v-if="gseOptions.length" class="sel" v-model="gse" @change="onGseChange">
          <option value="" disabled>选择 GSE…</option>
          <option v-for="g in gseOptions" :key="g" :value="g">{{ g }}</option>
        </select>
        <div v-else class="field">
          <input v-model="gse" placeholder="GSE63067" spellcheck="false"
                 @change="onGseChange" /><i class="line"></i>
        </div>
        <span v-if="gse" class="state on">已确认 {{ gse }}</span>
      </div>

      <!-- 预检展示 -->
      <template v-if="tab !== 'local' && gse && (probeBusy || probe)">
        <p v-if="probeBusy" class="hint-line">预检中(测速 + 可用性判定)…</p>
        <template v-else>
          <p v-if="probeR" class="hint-line">
            矩阵 {{ fmtSize(probeR.size_bytes) }} · 实测 {{ (probeR.speed_bps / 1024).toFixed(0) }} KB/s ·
            预计 {{ fmtEta(Math.min(probeR.size_bytes / Math.max(probeR.speed_bps, 1), 5400)) }}
          </p>
          <p v-if="probeR && (!probeR.reliable || probeR.size_bytes / Math.max(probeR.speed_bps, 1) > 300)"
             class="stale">带宽不稳或预计超 5 分钟 —— 慢速时数据仍在流动,请耐心等待。</p>
          <p v-if="probeErr" class="hint-line">预检失败(不影响正式运行):{{ probeErr }}</p>
          <template v-if="probeAssess">
            <p v-if="probeAssess.matrix_empty && probeAssess.counts_available" class="hint-line good">
              series matrix 为元数据空表 —— 检测到 counts 补充文件,运行时自动回退 DESeq2,可用。
            </p>
            <p v-else-if="probeAssess.matrix_empty" class="stale">
              空矩阵且无 counts 补充 —— 该 GSE 确定性不可用,配置与运行已锁定。
            </p>
            <p v-else-if="!probeAssess.symbol_level && probeAssess.annot === 'missing'"
               class="stale">探针级矩阵且平台注释缺失 —— 确定性不可用,已锁定。</p>
          </template>
          <p v-else-if="probe?.annot_missing" class="stale">
            平台 {{ probe?.gpl }} 注释文件不可达 —— 大概率不可用,总览后做权威判定。
          </p>
        </template>
      </template>
    </section>

    <!-- ============ ② 配置区 ============ -->
    <section class="sec">
      <h3 class="sec-title">② 配置</h3>
      <p v-if="unusable" class="stale">数据集已被判定不可用,配置与运行已锁定。</p>
      <div class="set-row">
        <label class="set-label">case 正则</label>
        <div class="field wide"><input v-model="casePat" :disabled="unusable"
                                       placeholder="steatosis|non-alcoholic" spellcheck="false" /><i class="line"></i></div>
      </div>
      <div class="set-row">
        <label class="set-label">control 正则</label>
        <div class="field wide"><input v-model="ctrlPat" :disabled="unusable"
                                       placeholder="healthy" spellcheck="false" /><i class="line"></i></div>
        <button class="bar-btn" :disabled="!gse || overviewJob.running.value || unusable"
                @click="loadOverview">
          {{ overviewJob.running.value ? '总览生成中…' : '生成分组总览(含下载)' }}
        </button>
      </div>
      <div v-if="overviewJob.running.value" class="progress">
        <div class="progress-track"><div class="progress-fill indeterminate"></div></div>
        <p class="progress-text">{{ overviewJob.text.value || '下载数据集中…' }}</p>
      </div>
      <p v-if="overviewJob.error.value" class="stale">总览失败:{{ overviewJob.error.value }}</p>

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
          数据集问题:{{ overview.issues.join('; ') }}
        </p>
        <template v-if="overview.counts_rows?.length">
          <p class="hint-line good">counts 数据集(共 {{ overview.counts_rows.length }} 列)—— 正则的匹配对象是下列样本列名:</p>
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
            矩阵 {{ overview.assess.shape[0] }}×{{ overview.assess.shape[1] }} ·
            ID 级别:{{ overview.assess.symbol_level ? '基因符号' : '探针级' }}
            {{ overview.assess.symbol_level ? '' : (overview.assess.annot === 'cached' || overview.assess.annot === 'ok'
              ? `· 注释可用(${overview.assess.gpl})`
              : overview.assess.annot === 'missing' ? `· 注释缺失(${overview.assess.gpl}) —— 不可用`
                : `· 注释未知(${overview.assess.gpl})`) }}
          </p>
          <p v-if="overview.assess.group_fields?.length" class="hint-line">
            疑似分组字段:{{ overview.assess.group_fields.slice(0, 3).join('; ') }}
          </p>
          <p v-if="overview.assess.n_samples < 12" class="stale">
            样本量偏低:{{ overview.assess.n_samples }}(建议每组 ≥ 6)
          </p>
        </template>
      </template>

      <!-- 排序 CSV -->
      <div class="set-row" style="margin-top:16px">
        <label class="set-label">平台排序表</label>
        <div class="field"><input v-model="drugFilter" spellcheck="false"
                                  placeholder="按药物过滤(可留空)" @change="loadRanked" /><i class="line"></i></div>
        <span class="ctl-meta">{{ rankedShown.length }} 个候选 · 勾选 {{ rankedSel.size }}</span>
      </div>
      <div v-if="rankedShown.length" class="filters">
        <label v-for="f in rankedShown" :key="f.label" class="chk">
          <input type="checkbox" :checked="rankedSel.has(f.label)" @change="toggleRanked(f.label)"
                 :disabled="unusable" /> {{ f.label }}
        </label>
      </div>
      <p v-else-if="disease" class="hint-line">该疾病下没有可用的排序表。</p>
      <p v-else class="hint-line">填写疾病名称后自动发现平台排序表(本机分析产物)。</p>

      <div class="opt-grid">
        <label class="opt">
          <span class="opt-label">shortlist Top-K</span>
          <Stepper v-model="topK" :min="5" :max="100" />
        </label>
        <label class="opt">
          <span class="opt-label">DEG 预过滤 Top-N</span>
          <Stepper v-model="prefilterN" :min="50" :max="5000" :step="50" />
        </label>
        <label class="opt">
          <span class="opt-label">Bootstrap 轮数</span>
          <Stepper v-model="bootstrapN" :min="5" :max="200" :step="5" />
        </label>
      </div>
      <p class="hint-line">默认 500 基因 × 30 轮:ML 段本机实测约 2 分钟;SVM-RFE 对特征数近平方级,调大会显著变慢。</p>
    </section>

    <!-- ============ ③ 运行区 ============ -->
    <section class="sec">
      <h3 class="sec-title">③ 运行</h3>
      <p v-if="unusable" class="stale">数据集不可用,运行已锁定。</p>
      <p v-else-if="sourceConflict" class="stale">GSE 与本地矩阵二选一 —— 请清掉其中一路。</p>
      <p v-else-if="localPartial" class="stale">本地矩阵:表达谱与标签两路必须同时填写。</p>
      <button class="build" :disabled="!canRun || runJob.running.value || overviewJob.running.value"
              @click="run">
        <svg viewBox="0 0 24 14" width="26" height="15">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <span class="build-hint">开始筛选</span>

      <div v-if="runJob.running.value" class="prog">
        <div class="progress">
          <div class="progress-track">
            <div class="progress-fill"
                 :style="{ width: `${((runJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
          </div>
        </div>
        <ul class="stages">
          <li v-for="k in stageOrder" :key="k" class="stage" :class="stageStates[k]">
            <i class="stage-icon">{{ STATE_ICON[stageStates[k]] }}</i>{{ STAGE_ZH[k] }}
          </li>
        </ul>
        <p v-if="runJob.text.value" class="hint-line">{{ runJob.text.value }}</p>
      </div>
      <p v-if="runJob.error.value" class="stale">筛选失败:{{ runJob.error.value }}</p>

      <div v-if="fallback" class="set-row" style="margin-top:14px">
        <span class="stale">⚠ GEO 下载失败,已回退本地缓存(缓存日期 {{ fallback.date }})—— 建议强制重下:</span>
        <button class="bar-btn" :disabled="forceJob.running.value" @click="forceRedownload">
          {{ forceJob.running.value ? '重下中…' : '强制重下' }}
        </button>
      </div>
      <p v-if="forceJob.error.value" class="stale">重下失败:{{ forceJob.error.value }}</p>
    </section>

    <!-- ============ ④ 结果区 ============ -->
    <section class="sec" v-if="resultData">
      <h3 class="sec-title">④ 结果({{ resultData.gse }} · {{ resultData.elapsed }}s)</h3>
      <div class="metrics">
        <div class="metric"><b>{{ meta?.n_samples }}</b><span>样本</span></div>
        <div class="metric"><b>{{ meta?.n_genes }}</b><span>基因</span></div>
        <div class="metric"><b>{{ meta?.cv_auc ?? 'N/A' }}</b><span>CV-AUC</span></div>
        <div class="metric"><b>{{ meta?.algo_spearman_mean }}</b><span>Spearman 均值</span></div>
      </div>

      <template v-if="resultData.roc_png_rel">
        <p class="hint-line" style="margin-top:16px">交叉验证 ROC(产物 PNG,论文可直接取用):</p>
        <img class="fig-img" :src="artifactsUrl(resultData.roc_png_rel)" alt="ROC" />
      </template>

      <h4 class="sub-sec">共识表(Top {{ Math.min(20, consensusRows.length) }} / 共 {{ consensusTotal }},全量见下方 CSV 下载)</h4>
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

      <h4 class="sub-sec">漏斗</h4>
      <div class="tbl-wrap" style="max-height:240px">
        <table class="tbl">
          <thead><tr><th v-for="k in Object.keys(funnelRows[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
          <tbody><tr v-for="(r, i) in funnelRows" :key="i"><td v-for="(v, k) in r" :key="k">{{ v }}</td></tr></tbody>
        </table>
      </div>

      <h4 class="sub-sec">Shortlist({{ shortlistCount }} 个)</h4>
      <div class="tbl-wrap" style="max-height:320px">
        <table class="tbl">
          <thead><tr><th v-for="k in Object.keys(shortlistRows[0] ?? {})" :key="k">{{ k }}</th></tr></thead>
          <tbody><tr v-for="(r, i) in shortlistRows" :key="i"><td v-for="(v, k) in r" :key="k">{{ v }}</td></tr></tbody>
        </table>
      </div>

      <div class="bar">
        <a v-for="(rel, name) in artifactsMap" :key="name" class="bar-btn"
           :href="artifactsUrl(rel, true)" :download="name">{{ name }}</a>
      </div>
    </section>
    <section class="sec" v-else>
      <h3 class="sec-title">④ 结果</h3>
      <p class="hint-line">运行筛选后,这里展示共识表 / ROC / shortlist,产物同时落盘。</p>
    </section>

    <!-- ============ 历史结果(近 3 次) ============ -->
    <section class="sec">
      <h3 class="sec-title">历史结果(近 3 次)</h3>
      <p v-if="!history.length" class="hint-line">尚无历史结果;完成一次分析后,这里会显示最近 3 次。</p>
      <div v-else class="set-row">
        <select class="sel inline" v-model="histPick">
          <option v-for="r in history" :key="r.name" :value="r.name">{{ r.name }}</option>
        </select>
        <span v-if="histRow" class="ctl-meta">
          AUC {{ histRow.auc ?? 'N/A' }} · 样本 {{ histRow.n_samples ?? '?' }} ·
          shortlist {{ histRow.shortlist_count ?? '?' }} · v{{ histRow.ml_version ?? '?' }}
        </span>
        <span v-for="s in histRow?.stale ?? []" :key="s.type" class="stale">
          {{ s.type === 'version' ? `算法版本已更新(${s.old} → ${s.new})`
            : s.type === 'ranked_changed' ? `排序表已变更:${s.name}`
              : `排序表缺失:${s.name}` }}
        </span>
      </div>
    </section>

    <p class="hint-line">
      治理与批量操作(删除 / 打包 / keep-N)在
      <a class="bar-btn" href="#/ml-results">结果管理页</a>
    </p>
  </div>
</template>

<style scoped>
.set-row {
  display: flex;
  align-items: baseline;
  gap: 18px;
  padding: 10px 0;
  flex-wrap: wrap;
}
.set-label {
  flex: none;
  width: 120px;
  font-size: 11.5px;
  letter-spacing: 0.1em;
  color: #8a8a88;
}
.state.on {
  color: #1a1a1a;
  font-weight: 500;
}
.sel.inline {
  width: auto;
  min-width: 200px;
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
.metrics {
  display: flex;
  gap: 56px;
  margin: 8px 0 6px;
}
.metric {
  display: flex;
  flex-direction: column;
  gap: 2px;
}
.metric b {
  font-size: 28px;
  font-weight: 300;
  font-variant-numeric: tabular-nums;
  color: #1a1a1a;
}
.metric span {
  font-size: 10.5px;
  letter-spacing: 0.2em;
  color: #bcbcbc;
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
.stale {
  font-size: 12px;
  letter-spacing: 0.05em;
  color: #b3261e;
}
.hint-line.good {
  color: #1a1a1a;
}
</style>
