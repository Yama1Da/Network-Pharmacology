<script setup lang="ts">
/** 主页(App.py 职责的 TS 化):四库检索落地页,1:1 复刻 Front-end-screen-options 极简主义样本。
 * 箭头点击 → /api/targets/retrieve job(四库并发检索,写入来源与 merged 缓存),
 * 完成后展示两侧靶点数并可预填跳转分析页。 */
import { computed, onBeforeUnmount, onMounted, reactive, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from '../api'
import { setPrefill } from '../store'
import { useJob } from '../composables/useJob'
import { theme } from '../utils/theme'

const { t } = useI18n()
const drug = ref('metformin')
const disease = ref('Type 2 diabetes mellitus')

/* 检索结果三张靶点表:药物 / 疾病 / 交集,页签切换(merged 表的行,来源标记 + 各库分数) */
type TargetRow = Record<string, string | number | null>
type ListTab = 'drug' | 'disease' | 'common'
const listTab = ref<ListTab>('common')

const job = useJob((result) => {
  done.value = {
    drug: Number(result.drug_count ?? 0),
    disease: Number(result.disease_count ?? 0),
    common_count: Number(result.common_count ?? 0),
    common_genes: (result.common_genes ?? []) as string[],
    failed: (result.failed ?? []) as string[],
    rows: {
      drug: (result.drug_rows ?? []) as TargetRow[],
      disease: (result.disease_rows ?? []) as TargetRow[],
      common: (result.common_rows ?? []) as TargetRow[]
    }
  }
  listTab.value = 'common'
})

const done = ref<{
  drug: number
  disease: number
  common_count: number
  common_genes: string[]
  failed: string[]
  rows: Record<ListTab, TargetRow[]>
} | null>(null)

type ColKind = 'sym' | 'flag' | 'score'
interface ListCol { key: string; label: string; kind: ColKind }

const listCols = computed<Record<ListTab, ListCol[]>>(() => {
  const drug: ListCol[] = [
    { key: 'approvedSymbol', label: t('web.home.col_gene'), kind: 'sym' },
    { key: 'source_chmbl', label: 'ChEMBL', kind: 'flag' },
    { key: 'source_dgidb', label: 'DGIdb', kind: 'flag' },
    { key: 'score_pchembl', label: 'pChEMBL', kind: 'score' },
    { key: 'score_dgidb', label: t('web.home.col_dgidb_score'), kind: 'score' }
  ]
  const disease: ListCol[] = [
    { key: 'approvedSymbol', label: t('web.home.col_gene'), kind: 'sym' },
    { key: 'source_ot', label: 'Open Targets', kind: 'flag' },
    { key: 'source_disgenet', label: 'DisGeNET', kind: 'flag' },
    { key: 'score_ot', label: t('web.home.col_ot_score'), kind: 'score' },
    { key: 'score_disgenet', label: t('web.home.col_disgenet_score'), kind: 'score' }
  ]
  return { drug, disease, common: [...drug, ...disease.slice(1)] }
})

const activeRows = computed<TargetRow[]>(() => done.value?.rows[listTab.value] ?? [])

function fmtCell(v: string | number | null | undefined, kind: ColKind): string {
  if (kind === 'flag') return v ? '✓' : '–'
  if (kind === 'score') return v == null ? '—' : Number(v).toFixed(3)
  return String(v ?? '')
}

function retrieve(): void {
  if (!drug.value.trim() || !disease.value.trim()) {
    return
  }
  done.value = null
  void job.run(api.retrieveTargets(drug.value.trim(), disease.value.trim()),
    { doneText: t('web.home.retrieve_done') })
}

function goAnalysis(): void {
  setPrefill(drug.value.trim(), disease.value.trim())
  location.hash = '#/analysis'
}

/* 页脚统计数字滚动(样本同款):数值为 /api/stats 的本机缓存实数,非装饰值 */
type StatKey = 'databases' | 'compounds' | 'targets' | 'pathways'
const STAT_KEYS: StatKey[] = ['databases', 'compounds', 'targets', 'pathways']
const stats = ref<Record<StatKey, number>>(
  { databases: 0, compounds: 0, targets: 0, pathways: 0 })
const shown = reactive<Record<StatKey, number>>(
  { databases: 0, compounds: 0, targets: 0, pathways: 0 })
const fmt = (n: number): string => n.toLocaleString('en-US')

let raf = 0
let started = false
let introDone = false
let statsDone = false
function maybeStart(): void {
  if (started || !introDone || !statsDone) return
  started = true
  const t0 = performance.now()
  const dur = 1600
  const step = (t: number): void => {
    const p = Math.min(1, (t - t0) / dur)
    const e = 1 - Math.pow(1 - p, 3)
    for (const k of STAT_KEYS) shown[k] = Math.round(stats.value[k] * e)
    if (p < 1) raf = window.requestAnimationFrame(step)
  }
  raf = window.requestAnimationFrame(step)
}
onMounted(() => {
  window.setTimeout(() => { introDone = true; maybeStart() }, 1750)  // 与 App 入场幕布对齐
  void api.homeStats()
    .then((s) => { stats.value = s })
    .catch(() => { /* 读取失败保持 0:宁可真实为 0,不放假数 */ })
    .finally(() => { statsDone = true; maybeStart() })
})
onBeforeUnmount(() => window.cancelAnimationFrame(raf))
</script>

<template>
  <p class="kicker s1">{{ t('web.home.kicker') }}</p>
  <h1 class="title s2">{{ t('web.home.title') }}</h1>
  <p class="sub s3">{{ t('web.home.sub') }}</p>

  <div class="home-form s4">
    <label class="field">
      <input v-model="drug" spellcheck="false" :placeholder="t('web.common.ph_drug')" />
      <i class="line"></i>
    </label>
    <span class="cross">×</span>
    <label class="field">
      <input v-model="disease" spellcheck="false" :placeholder="t('web.common.ph_disease')" />
      <i class="line"></i>
    </label>
    <button class="go" :disabled="job.running.value" :title="t('web.home.retrieve_title')"
            :aria-label="t('web.home.retrieve_title')"
            @click="retrieve">
      <svg v-if="!job.running.value" viewBox="0 0 24 14" width="22" height="13">
        <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
              stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
      </svg>
      <span v-else class="spin"></span>
    </button>
  </div>

  <Transition name="note">
    <div v-if="job.running.value" class="retrieve">
      <div class="progress">
        <div class="progress-track">
          <div class="progress-fill"
               :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
        </div>
        <p class="progress-text">{{ job.text.value || t('web.home.retrieving') }}</p>
      </div>
    </div>
  </Transition>

  <p v-if="job.error.value" class="retrieve failed-msg">
    {{ t('web.home.retrieve_failed', { msg: job.error.value }) }}
  </p>

  <Transition name="note">
    <div v-if="done" class="retrieve done">
      <p class="done-line">
        {{ t('web.home.done_drug') }} <b>{{ done.drug }}</b> ·
        {{ t('web.home.done_disease') }} <b>{{ done.disease }}</b> ·
        {{ t('web.home.done_common') }} <b>{{ done.common_count }}</b>
        <span v-if="done.failed.length" class="failed">{{ t('web.home.failed_sources', { names: done.failed.join(', ') }) }}</span>
      </p>
      <img v-if="job.jobId.value" class="venn" :src="`/api/targets/venn?job_id=${job.jobId.value}&theme=${theme}`"
           :alt="t('web.home.venn_alt')" />
      <div class="tl">
        <div class="tl-tabs">
          <button class="lv" :class="{ on: listTab === 'drug' }" @click="listTab = 'drug'">
            {{ t('web.home.tab_drug') }}&nbsp;{{ done.drug }}</button>
          <span class="lv-sep">/</span>
          <button class="lv" :class="{ on: listTab === 'disease' }" @click="listTab = 'disease'">
            {{ t('web.home.tab_disease') }}&nbsp;{{ done.disease }}</button>
          <span class="lv-sep">/</span>
          <button class="lv" :class="{ on: listTab === 'common' }" @click="listTab = 'common'">
            {{ t('web.home.tab_common') }}&nbsp;{{ done.common_count }}</button>
        </div>
        <div class="tl-scroll">
          <table v-if="activeRows.length" class="tbl">
            <thead>
              <tr>
                <th v-for="c in listCols[listTab]" :key="c.key">{{ c.label }}</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(r, i) in activeRows" :key="i">
                <td v-for="c in listCols[listTab]" :key="c.key"
                    :class="c.kind === 'sym' ? 'sym' : c.kind === 'flag' ? { flag: true, on: !!r[c.key] } : 'score'">
                  {{ fmtCell(r[c.key], c.kind) }}</td>
              </tr>
            </tbody>
          </table>
          <p v-else class="tl-empty">{{ t('web.home.list_empty') }}</p>
        </div>
      </div>
      <button class="go enter" :title="t('web.home.enter_title')" @click="goAnalysis">
        <svg viewBox="0 0 24 14" width="22" height="13">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <p class="hint-line">{{ t('web.home.ready_hint') }}</p>
    </div>
  </Transition>

  <div class="foot-stats s5">
    <span><b>{{ fmt(shown.databases) }}</b> {{ t('web.home.stat_dbs') }} ·
      <b>{{ fmt(shown.compounds) }}</b> {{ t('web.home.stat_compounds') }} ·
      <b>{{ fmt(shown.targets) }}</b> {{ t('web.home.stat_targets') }} ·
      <b>{{ fmt(shown.pathways) }}</b> {{ t('web.home.stat_pathways') }}</span>
  </div>
</template>

<style scoped>
.home-form {
  margin-top: 56px;
  display: flex;
  align-items: flex-end;
  gap: 26px;
  flex-wrap: wrap;
  justify-content: center;
}
.field {
  position: relative;
  width: clamp(180px, 17vw, 250px);
}
.field input {
  width: 100%;
  box-sizing: border-box;
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid var(--c-line);
  padding: 10px 2px;
  font-size: 15.5px;
  font-family: inherit;
  color: var(--c-ink);
  letter-spacing: 0.02em;
  border-radius: 0;
}
.field input::placeholder {
  color: var(--c-faint);
  letter-spacing: 0.08em;
}
.line {
  position: absolute;
  left: 0;
  right: 0;
  bottom: 0;
  height: 1px;
  background: var(--c-accent);
  transform: scaleX(0);
  transition: transform 0.45s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.field:focus-within .line {
  transform: scaleX(1);
}
.cross {
  padding-bottom: 12px;
  color: var(--c-faint2);
  font-size: 15px;
}
.go {
  width: 46px;
  height: 46px;
  margin-bottom: 1px;
  border: none;
  border-radius: 50%;
  background: var(--c-accent);
  color: var(--c-on-accent);
  display: grid;
  place-items: center;
  cursor: pointer;
  font-family: inherit;
  transition: transform 0.3s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.go:hover:not(:disabled) {
  transform: scale(1.06);
}
.go:disabled {
  opacity: 0.55;
  cursor: default;
}
.spin {
  width: 14px;
  height: 14px;
  border: 1.6px solid color-mix(in srgb, var(--c-on-accent) 35%, transparent);
  border-top-color: var(--c-on-accent);
  border-radius: 50%;
  animation: rot 0.8s linear infinite;
}
@keyframes rot {
  to {
    transform: rotate(360deg);
  }
}
.retrieve {
  /* 封顶用 % 不用 vw:vw 不随界面缩放(CSS zoom)收缩,放大后会出现横向滚动条 */
  width: min(560px, 90%);
  margin-top: 34px;
}
/* 完成态装得下三页签靶点表,容器独立加宽;.tl 直接 100% 不再反向 breakout */
.retrieve.done {
  width: min(880px, 90%);
}
.failed-msg {
  margin-left: auto;
  margin-right: auto;
  font-size: 12px;
  letter-spacing: 0.05em;
  color: var(--c-danger);
  text-align: center;
}
.done {
  text-align: center;
}
.done-line {
  margin: 0 0 14px;
  font-size: 13.5px;
  letter-spacing: 0.06em;
  color: var(--c-ink);
}
.done-line b {
  font-weight: 600;
  font-variant-numeric: tabular-nums;
}
.failed {
  font-size: 11px;
  color: var(--c-danger);
}
.enter {
  margin: 0 auto;
}
.venn {
  max-width: min(520px, 90%);
  margin: 6px auto 0;
  display: block;
  background: var(--c-panel);
  border: 1px solid var(--c-line-soft);
}
.tl {
  width: 100%;
  margin-top: 20px;
}
.tl-tabs {
  display: flex;
  justify-content: center;
  gap: 14px;
  margin-bottom: 12px;
}
.tl-scroll {
  max-height: 280px;
  overflow: auto;
  border: 1px solid var(--c-line-soft);
  padding: 0 14px;
}
.tl-scroll .tbl th,
.tl-scroll .tbl td {
  padding-right: 12px;
}
.tl-scroll td.sym {
  font-weight: 500;
  color: var(--c-ink);
}
.tl-scroll td.flag {
  color: var(--c-flag-off);
}
.tl-scroll td.flag.on {
  color: var(--c-ink);
}
.tl-scroll td.score {
  color: var(--c-mid2);
}
.tl-empty {
  margin: 20px 0;
  font-size: 11.5px;
  letter-spacing: 0.05em;
  color: var(--c-ink3);
}
.foot-stats {
  margin-top: 48px;
  font-size: 11px;
  letter-spacing: 0.1em;
  color: var(--c-ink3);
}
.foot-stats b {
  font-weight: 500;
  color: var(--c-num2);
  font-variant-numeric: tabular-nums;
}
.hint-line {
  margin: 10px 0 0;
  font-size: 11px;
  letter-spacing: 0.06em;
  color: var(--c-ink4);
}
</style>
