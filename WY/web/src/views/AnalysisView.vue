<script setup lang="ts">
/** PPI 分析(Data_Network 的 TS 迁移,阶段 3a)。
 * 产物直出:分析图表由后端落盘、经 /api/artifacts 伺服;交互区重跑(模块/GSEA)
 * 无状态化——从落盘 nodes/edges CSV 重建图。行为与 Streamlit 版一致。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import {
  api, artifactsUrl, type PpiAnalyzeBody, type PpiResult, type RwrRow, type Row
} from '../api'
import { takePrefill } from '../store'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import Stepper from '../components/Stepper.vue'
import { downloadBlob } from '../utils/download'
import { fmt } from '../utils/format'

const { t } = useI18n()

/* ---------------- 表单 ---------------- */

const mode = ref<'names' | 'genes'>('names')
const drug = ref('metformin')
const disease = ref('Type 2 diabetes mellitus')
const genesText = ref('')
const score = ref(0.45)
const networkType = ref('full')
const rwrExpand = ref(0)
const hubTopN = ref(10)
const hubMethod = ref('degree')

const snapshot = computed<(string | number)[]>(() => [
  mode.value, drug.value.trim(), disease.value.trim(), genesText.value.trim(),
  score.value, networkType.value, rwrExpand.value, hubTopN.value, hubMethod.value
])

/* ---------------- 分析任务 ---------------- */

const ppi = ref<PpiResult | null>(null)
const job = useJob((result) => {
  ppi.value = result as unknown as PpiResult
})

const stages = computed<string[]>(() => {
  const base = ['targets', 'string', 'graph', 'enrichment']
  const extra = mode.value === 'names' ? ['gsea', 'network'] : []
  return [...base, ...extra, 'modules', 'save']
})

const stageLabels = computed<Record<string, string>>(() => ({
  targets: t('web.analysis.stage_targets'),
  string: t('web.analysis.stage_string'),
  graph: t('web.analysis.stage_graph'),
  enrichment: t('web.analysis.stage_enrichment'),
  gsea: t('web.analysis.stage_gsea'),
  network: t('web.analysis.stage_network'),
  modules: t('web.analysis.stage_modules'),
  save: t('web.analysis.stage_save')
}))
const STATE_ICON: Record<string, string> = { pending: '○', running: '◐', done: '●' }

const stepStates = computed<Record<string, string>>(() => {
  const map: Record<string, string> = {}
  for (const k of stages.value) map[k] = 'pending'
  const p = job.extra.value as { step?: string; state?: string; steps?: Record<string, string> }
  // 后端富集/GSEA/邻近度+RWR/模块检测并发执行,整体上报 steps 映射 → 多阶段同时点亮;
  // steps 缺失时(旧任务/其他作业)退回单 step/state 语义
  if (p?.steps) {
    for (const [k, st] of Object.entries(p.steps)) if (k in map) map[k] = st
  } else if (p?.step && p?.state && p.step in map) {
    map[p.step] = p.state
  }
  return map
})

function analyze(): void {
  if (mode.value === 'names' && (!drug.value.trim() || !disease.value.trim())) {
    flash(t('web.analysis.need_names'))
    return
  }
  const body: PpiAnalyzeBody = {
    mode: mode.value, drug: drug.value, disease: disease.value, genes_text: genesText.value,
    score: score.value, network_type: networkType.value, rwr_expand: rwrExpand.value,
    hub_top_n: hubTopN.value, hub_method: hubMethod.value
  }
  void job.run(api.analyzePpi(body), { doneText: t('web.analysis.report_done') })
}

/* ---------------- 结果 ---------------- */

const paramsStale = computed(() => {
  if (!ppi.value) return false
  const a = snapshot.value
  const b = ppi.value.params
  return a.length !== b.length || a.some((v, i) => String(v) !== String(b[i]))
})

const { notice, flash } = useNotice(3200)

const view = ref('')
const availableViews = computed(() => ppi.value?.views ?? [])
const effectiveView = computed(() => {
  const v = view.value
  if (v && availableViews.value.includes(v)) return v
  return availableViews.value.includes('cytoscape') ? 'cytoscape' : 'interactive'
})
function setView(v: string): void {
  view.value = v
}
const viewLabels = computed<Record<string, string>>(() => ({
  interactive: t('web.analysis.view_interactive'),
  static: t('web.analysis.view_static'),
  string: t('web.analysis.view_string'),
  cytoscape: t('web.analysis.view_cytoscape')
}))

const artifacts = computed(() => ppi.value?.artifacts ?? {})
function art(key: string): string | null {
  return artifacts.value[key] ?? null
}

/* 邻近度(PpiResult.proximity 已是 ProximityResult,无需再转 any) */
const prox = computed(() => ppi.value?.proximity ?? null)

/* RWR */
const rwr = computed(() => ppi.value?.rwr ?? null)
const rwrSide = ref<'disease_seeds' | 'drug_seeds'>('disease_seeds')
const rwrRes = computed(() => rwr.value?.[rwrSide.value] ?? null)
const rwrTop = computed<(RwrRow & { is_seed: boolean })[]>(() => {
  const res = rwrRes.value
  if (!res?.ok) return []
  const seeds = (res.seeds ?? []).map((r) => ({ ...r, is_seed: true }))
  const non = (res.top_nonseed ?? []).map((r) => ({ ...r, is_seed: false }))
  return [...seeds, ...non].sort((a, b) => b.score - a.score)
})
const rwrNonSeed = computed(() => rwrTop.value.filter((r) => !r.is_seed))
const rwrSeeds = computed(() => rwrTop.value.filter((r) => r.is_seed))

/* 模块 / GSEA 重跑 */
const modsJob = useJob((result) => {
  if (ppi.value) {
    ppi.value = { ...ppi.value, modules: {
      summary: (result.summary ?? []) as Row[],
      failures: (result.failures ?? []) as unknown[],
      no_sig: (result.no_sig ?? []) as number[],
      hidden_single: Number(result.hidden_single ?? 0),
      enr_rows: (result.enr_rows ?? []) as Row[],
      enr_total: Number(result.enr_total ?? 0),
      enr_modules: Number(result.enr_modules ?? 0)
    } }
  }
  flash(t('web.analysis.modules_done_flash'))
})
const gseaJob = useJob((result) => {
  if (ppi.value && gseaSide.value) {
    const g = { ...ppi.value.gsea }
    g[gseaSide.value] = {
      reason: String(result.reason ?? ''),
      stats: (result.stats ?? {}) as Row,
      rows: (result.rows ?? []) as Row[]
    }
    ppi.value = { ...ppi.value, gsea: g }
  }
  flash(t('web.analysis.gsea_done_flash'))
})
const gseaSide = ref<'drug' | 'disease'>('drug')
const gseaRes = computed(() => ppi.value?.gsea?.[gseaSide.value] ?? null)
const gseaDone = computed(() => {
  const g = gseaRes.value
  return !!g && !paramsStale.value &&
    !['request_failed', 'circuit_open', 'normalize_empty'].includes(g.reason)
})

/* 下载清单 */
const downloadItems = computed<{ label: string; rel: string }[]>(() => {
  const a = artifacts.value
  const out: { label: string; rel: string }[] = []
  for (const key of ['nodes', 'edges', 'html', 'graphml', 'enrichment_csv',
    'network_png', 'network_svg', 'string_png', 'string_svg', 'cyto_png', 'cyto_svg',
    'core_png', 'core_svg', 'dtp_png', 'dtp_svg',
    'enrichment_bar_png', 'enrichment_bar_svg', 'enrichment_dot_png', 'enrichment_dot_svg',
    'gsea_drug_png', 'gsea_drug_svg', 'gsea_disease_png', 'gsea_disease_svg',
    'gsea_drug_csv', 'gsea_disease_csv', 'rwr_disease_csv', 'rwr_drug_csv',
    'proximity_csv', 'query_genes_csv', 'modules_csv', 'module_enrichment_csv']) {
    if (a[key]) out.push({ label: key, rel: a[key] })
  }
  return out
})

const dirName = computed(() =>
  (ppi.value?.directory_rel ?? '').split('/').pop() ?? '')

async function zipWholeDirectory(): Promise<void> {
  if (!dirName.value) return
  try {
    downloadBlob(await api.zip('downloads', [dirName.value]), `${dirName.value}.zip`)
  } catch {
    flash(t('web.analysis.dl_zip_fail'))
  }
}

onMounted(() => {
  // 主页四库检索完成后的预填跳转
  const pf = takePrefill()
  if (pf) {
    mode.value = 'names'
    drug.value = pf.drug
    disease.value = pf.disease
  }
})
</script>

<template>
  <p class="kicker s1">{{ t('web.analysis.kicker') }}</p>
  <h1 class="title s2">{{ t('web.analysis.title') }}</h1>
  <p class="sub s3">{{ t('web.analysis.sub') }}</p>

  <!-- ================= 表单 ================= -->
  <div class="panel wide s4">
    <div class="xlinks">
      <a class="bar-btn" href="#/data">{{ t('web.nav.data_link') }}</a>
    </div>
    <div class="mode-row">
      <span class="ctl-label">{{ t('web.analysis.mode_label') }}</span>
      <span class="lv-group">
        <button class="lv" :class="{ on: mode === 'names' }" @click="mode = 'names'">{{ t('web.analysis.mode_names') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: mode === 'genes' }" @click="mode = 'genes'">{{ t('web.analysis.mode_genes') }}</button>
      </span>
    </div>

    <div v-if="mode === 'names'" class="names-form">
      <label class="field">
        <input v-model="drug" spellcheck="false" :placeholder="t('web.common.ph_drug')" />
        <i class="line"></i>
      </label>
      <span class="cross">×</span>
      <label class="field">
        <input v-model="disease" spellcheck="false" :placeholder="t('web.common.ph_disease')" />
        <i class="line"></i>
      </label>
    </div>
    <label v-else class="field genes">
      <textarea v-model="genesText" rows="2" spellcheck="false"
                :placeholder="t('web.analysis.ph_genes')" />
      <i class="line"></i>
    </label>

    <div class="opt-grid">
      <label class="opt">
        <span class="opt-label">{{ t('web.analysis.opt_score', { v: score.toFixed(2) }) }}</span>
        <input v-model.number="score" type="range" min="0.15" max="0.9" step="0.05" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.analysis.opt_network') }}</span>
        <select v-model="networkType">
          <option value="full">{{ t('web.analysis.net_full') }}</option>
          <option value="physical">{{ t('web.analysis.net_physical') }}</option>
        </select>
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.analysis.opt_rwr') }}</span>
        <Stepper v-model="rwrExpand" :min="0" :max="100" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.analysis.opt_hub_n') }}</span>
        <Stepper v-model="hubTopN" :min="5" :max="20" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.analysis.opt_hub_method') }}</span>
        <select v-model="hubMethod">
          <option value="degree">{{ t('web.analysis.hub_degree') }}</option>
          <option value="mcc">{{ t('web.analysis.hub_mcc') }}</option>
        </select>
      </label>
    </div>

    <div class="build-row">
      <button class="build" :disabled="job.running.value || modsJob.running.value || gseaJob.running.value"
              :title="t('web.analysis.build_title')" :aria-label="t('web.analysis.build_title')" @click="analyze">
        <svg viewBox="0 0 24 14" width="26" height="15">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <span class="build-hint">{{ t('web.analysis.build_title') }}</span>
    </div>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <!-- 进度:阶段图标行 + 细线 -->
    <div v-if="job.running.value" class="prog">
      <div class="progress">
        <div class="progress-track">
          <div class="progress-fill"
               :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
        </div>
      </div>
      <ul class="stages">
        <li v-for="k in stages" :key="k" class="stage" :class="stepStates[k]">
          <i class="stage-icon">{{ STATE_ICON[stepStates[k]] }}</i>{{ stageLabels[k] }}
        </li>
      </ul>
    </div>
    <p v-if="job.error.value" class="notice err">{{ t('web.analysis.err', { msg: job.error.value }) }}</p>
  </div>

  <!-- ================= 结果 ================= -->
  <div v-if="ppi" class="panel wide result">
    <p v-if="paramsStale" class="stale">{{ t('web.analysis.params_stale') }}</p>

    <div class="metrics s5">
      <div class="metric"><b>{{ ppi.genes_count }}</b><span>{{ t('web.analysis.m_input') }}</span></div>
      <div class="metric"><b>{{ ppi.nodes }}</b><span>{{ t('web.analysis.m_nodes') }}</span></div>
      <div class="metric"><b>{{ ppi.edges }}</b><span>{{ t('web.analysis.m_edges') }}</span></div>
      <div class="metric"><b>{{ ppi.hubs }}</b><span>{{ t('web.analysis.m_hubs') }}</span></div>
    </div>

    <p v-if="ppi.expansion" class="hint-line">
      {{ t('web.analysis.rwr_expansion', { seeds: ppi.expansion.seeds, expanded: ppi.expansion.expanded,
                                            total: ppi.expansion.seeds + ppi.expansion.expanded,
                                            prefix: ppi.prefix }) }}
    </p>
    <p v-if="!ppi.genes_count" class="stale">{{ t('web.analysis.no_common') }}</p>
    <p v-if="!ppi.edges" class="stale">{{ t('web.analysis.no_edges') }}</p>

    <!-- Hub 表 -->
    <section v-if="ppi.hub_table.length" class="sec">
      <h3 class="sec-title">
        {{ ppi.hub_rule.startsWith('top')
          ? t('web.analysis.hub_top_title', { n: ppi.hub_table.length, method: ppi.hub_method_eff })
          : t('web.analysis.hub_title') }}
      </h3>
      <div class="tbl-wrap">
        <table class="tbl">
          <thead>
            <tr><th v-for="k in Object.keys(ppi.hub_table[0])" :key="k">{{ k }}</th></tr>
          </thead>
          <tbody>
            <tr v-for="(r, i) in ppi.hub_table" :key="i">
              <td v-for="(v, k) in r" :key="k">{{ v }}</td>
            </tr>
          </tbody>
        </table>
      </div>
    </section>

    <!-- 网络图:视图切换 -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.analysis.sec_network') }}</h3>
      <div class="lv-group">
        <button v-for="v in ppi.views" :key="v" class="lv"
                :class="{ on: effectiveView === v }" @click="setView(v)">{{ viewLabels[v] }}</button>
      </div>
      <p v-if="effectiveView === 'string'" class="hint-line">{{ t('web.analysis.view_string_hint') }}</p>
      <p v-if="effectiveView === 'cytoscape'" class="hint-line">{{ t('web.analysis.view_cyto_hint') }}</p>
      <iframe v-if="effectiveView === 'interactive'"
              class="fig-frame" :src="artifactsUrl(art('html') ?? '')" />
      <img v-else-if="effectiveView === 'static'"
           class="fig-img" :src="artifactsUrl(art('png') ?? '')" alt="static network" />
      <img v-else-if="effectiveView === 'string'"
           class="fig-img" :src="artifactsUrl(art('string_png') ?? '')" alt="string style" />
      <img v-else class="fig-img" :src="artifactsUrl(art('cyto_png') ?? '')" alt="cytoscape style" />
      <p class="hint-line">{{ t('web.analysis.views_saved_hint') }}</p>
    </section>

    <!-- 核心靶点 -->
    <section v-if="art('core_png')" class="sec">
      <h3 class="sec-title">{{ t('web.analysis.sec_core') }}</h3>
      <img class="fig-img" :src="artifactsUrl(art('core_png') ?? '')" alt="core targets" />
    </section>

    <!-- 邻近度 -->
    <section v-if="prox && ppi.params[0] === 'names'" class="sec">
      <h3 class="sec-title">{{ t('web.analysis.sec_prox') }}</h3>
      <p class="hint-line">
        {{ t('web.analysis.prox_setup', { n: prox.n_perm, seed: prox.seed ?? '-' }) }}
      </p>
      <template v-if="prox.error === 'not_built'">
        <p class="stale">{{ t('web.analysis.prox_not_built') }}</p>
      </template>
      <template v-else-if="prox.error &&
                           ['low_coverage_drug', 'low_coverage_disease', 'empty_after_filter'].includes(prox.error)">
        <p class="stale">
          {{ t('web.analysis.prox_low_cov', { din: prox.coverage?.drug_in_lcc ?? 0, dtotal: prox.coverage?.drug_total ?? 0,
                                               iin: prox.coverage?.disease_in_lcc ?? 0, itotal: prox.coverage?.disease_total ?? 0 }) }}
        </p>
      </template>
      <template v-else-if="prox.error">
        <p class="stale">{{ t('web.analysis.prox_fail', { err: prox.error }) }}</p>
      </template>
      <template v-else-if="prox.ok">
        <div class="metrics three">
          <div class="metric"><b>{{ fmt(prox.d_c) }}</b><span>{{ t('web.analysis.m_dc') }}</span></div>
          <div class="metric"><b>{{ fmt(prox.z) }}</b><span>{{ t('web.analysis.m_z') }}</span></div>
          <div class="metric"><b>{{ fmt(prox.p_empirical, 4) }}</b><span>{{ t('web.analysis.m_p') }}</span></div>
        </div>
        <p class="hint-line" :class="{ good: typeof prox.z === 'number' && prox.z < 0 }">
          {{ typeof prox.z === 'number'
            ? (prox.z < 0 ? t('web.analysis.prox_support') : t('web.analysis.prox_unsupported'))
            : t('web.analysis.prox_degenerate') }}
        </p>
        <p class="hint-line">
          {{ t('web.analysis.prox_lcc', { nodes: prox.lcc_nodes, din: prox.coverage?.drug_in_lcc ?? 0,
                                          dtotal: prox.coverage?.drug_total ?? 0,
                                          iin: prox.coverage?.disease_in_lcc ?? 0,
                                          itotal: prox.coverage?.disease_total ?? 0 }) }}
          {{ prox.cached ? t('web.analysis.cached_hit') : '' }}
        </p>
        <template v-if="prox.d_c_excl != null">
          <p class="hint-line">{{ t('web.analysis.prox_excl', { n: prox.n_excl ?? 0 }) }}</p>
          <div class="metrics three">
            <div class="metric"><b>{{ fmt(prox.d_c_excl) }}</b><span>{{ t('web.analysis.m_dc_excl') }}</span></div>
            <div class="metric"><b>{{ fmt(prox.z_excl) }}</b><span>{{ t('web.analysis.m_z_excl') }}</span></div>
            <div class="metric"><b>{{ fmt(prox.p_excl, 4) }}</b><span>{{ t('web.analysis.m_p_excl') }}</span></div>
          </div>
        </template>
      </template>
    </section>

    <!-- RWR -->
    <section v-if="rwr && ppi.params[0] === 'names'" class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">{{ t('web.analysis.sec_rwr') }}</summary>
        <p class="hint-line">{{ t('web.analysis.rwr_hint') }}</p>
        <div class="lv-group">
          <button class="lv" :class="{ on: rwrSide === 'disease_seeds' }"
                  @click="rwrSide = 'disease_seeds'">{{ t('web.analysis.rwr_side_disease') }}</button>
          <i class="lv-sep">/</i>
          <button class="lv" :class="{ on: rwrSide === 'drug_seeds' }"
                  @click="rwrSide = 'drug_seeds'">{{ t('web.analysis.rwr_side_drug') }}</button>
        </div>
        <template v-if="rwrRes">
          <p v-if="rwrRes.error === 'not_built'" class="stale">{{ t('web.analysis.rwr_not_built') }}</p>
          <p v-else-if="rwrRes.error" class="stale">{{ t('web.analysis.rwr_fail', { err: rwrRes.error }) }}</p>
          <template v-else-if="rwrRes.ok">
            <p class="hint-line">
              {{ t('web.analysis.rwr_meta', { n: rwrRes.n_iter,
                                              conv: rwrRes.converged ? t('web.analysis.converged') : t('web.analysis.not_converged'),
                                              in: rwrRes.coverage?.seed_in_lcc ?? 0, total: rwrRes.coverage?.seed_total ?? 0 }) }}
              {{ rwrRes.cached ? t('web.analysis.cached_hit') : '' }}
            </p>
            <p v-if="rwrRes.enrichment && rwrRes.enrichment.z != null" class="hint-line"
               :class="{ good: rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05 }">
              {{ t('web.analysis.rwr_enrich', {
                n: rwrRes.enrichment.n_annotated, z: fmt(rwrRes.enrichment.z),
                p: fmt(rwrRes.enrichment.p, 4),
                sig: rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05
                  ? t('web.analysis.sig_yes') : t('web.analysis.sig_no') }) }}
            </p>
            <div v-if="rwrNonSeed.length" class="tbl-wrap">
              <table class="tbl">
                <thead><tr><th>gene</th><th>score</th><th v-if="'annotated' in rwrNonSeed[0]">{{ t('web.analysis.th_annot') }}</th></tr></thead>
                <tbody>
                  <tr v-for="r in rwrNonSeed" :key="String(r.gene)">
                    <td>{{ r.gene }}</td><td>{{ fmt(r.score, 6) }}</td>
                    <td v-if="'annotated' in r">{{ r.annotated ? '✓' : '' }}</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p v-else class="hint-line">{{ t('web.analysis.rwr_empty') }}</p>
            <details v-if="rwrSeeds.length" class="fold">
              <summary class="hint-line">{{ t('web.analysis.seeds_fold', { n: rwrSeeds.length }) }}</summary>
              <div class="tbl-wrap">
                <table class="tbl">
                  <thead><tr><th>gene</th><th>score</th></tr></thead>
                  <tbody><tr v-for="r in rwrSeeds" :key="String(r.gene)"><td>{{ r.gene }}</td><td>{{ fmt(r.score, 6) }}</td></tr></tbody>
                </table>
              </div>
            </details>
          </template>
        </template>
      </details>
    </section>

    <!-- 富集 -->
    <section class="sec">
      <h3 class="sec-title">{{ t('web.analysis.sec_enrich') }}</h3>
      <p v-if="ppi.enrichment.reason === 'both_sources_down' || ppi.enrichment.reason === 'gprofiler_down_no_fallback'"
         class="stale">
        {{ t('web.analysis.enrich_down') }}
      </p>
      <p v-else-if="!ppi.enrichment.rows.length" class="hint-line">{{ t('web.analysis.enrich_empty') }}</p>
      <template v-else>
        <p class="hint-line">
          {{ t('web.analysis.enrich_meta', { src: ppi.enrichment.source, n: ppi.enrichment.rows.length }) }}
        </p>
        <div class="tbl-wrap">
          <table class="tbl">
            <thead><tr><th>{{ t('web.analysis.th_source') }}</th><th>{{ t('web.analysis.th_term_id') }}</th><th>{{ t('web.analysis.th_term_name') }}</th><th>{{ t('web.analysis.th_adj_p') }}</th><th>{{ t('web.analysis.th_intersection') }}</th><th>{{ t('web.analysis.th_query') }}</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in ppi.enrichment.rows" :key="i">
                <td>{{ r.source }}</td><td>{{ r.term_id }}</td><td>{{ r.term_name }}</td>
                <td>{{ fmt(r.adjusted_p_value, 4) }}</td><td>{{ r.intersection_size }}</td><td>{{ r.query_size }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <template v-if="art('enrichment_bar_png')">
          <p class="hint-line">{{ t('web.analysis.enrich_bar_hint') }}</p>
          <img class="fig-img" :src="artifactsUrl(art('enrichment_bar_png') ?? '')" alt="enrichment bar" />
        </template>
        <template v-if="art('enrichment_dot_png')">
          <p class="hint-line">{{ t('web.analysis.enrich_dot_hint') }}</p>
          <img class="fig-img" :src="artifactsUrl(art('enrichment_dot_png') ?? '')" alt="enrichment dot" />
        </template>
        <template v-if="art('dtp_png')">
          <p class="hint-line">{{ t('web.analysis.dtp_hint') }}</p>
          <img class="fig-img" :src="artifactsUrl(art('dtp_png') ?? '')" alt="drug target pathway" />
        </template>
      </template>
    </section>

    <!-- GSEA -->
    <section v-if="ppi.params[0] === 'names'" class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">{{ t('web.analysis.sec_gsea') }}</summary>
        <div class="lv-group">
          <button class="lv" :class="{ on: gseaSide === 'drug' }" @click="gseaSide = 'drug'">{{ t('web.analysis.gsea_side_drug') }}</button>
          <i class="lv-sep">/</i>
          <button class="lv" :class="{ on: gseaSide === 'disease' }" @click="gseaSide = 'disease'">{{ t('web.analysis.gsea_side_disease') }}</button>
        </div>
        <button class="bar-btn" :disabled="gseaDone || paramsStale || gseaJob.running.value"
                @click="gseaSide && gseaJob.run(api.rerunGsea(gseaSide, drug, disease, ppi.prefix), { doneText: t('web.analysis.gsea_done_text') })">
          {{ gseaJob.running.value ? t('web.analysis.rerunning') : (gseaDone ? t('web.analysis.gsea_has_result') : t('web.analysis.gsea_rerun')) }}
        </button>
        <p v-if="gseaJob.running.value" class="hint-line">{{ gseaJob.text.value || t('web.analysis.rerunning') }}</p>
        <p v-if="gseaJob.error.value" class="notice err">{{ t('web.analysis.gsea_done_text') }}:{{ gseaJob.error.value }}</p>
        <template v-if="gseaRes">
          <p v-if="gseaRes.reason === 'insufficient_input'" class="stale">{{ t('web.analysis.gsea_insufficient') }}</p>
          <p v-else-if="gseaRes.reason === 'request_failed'" class="stale">{{ t('web.analysis.gsea_req_failed') }}</p>
          <p v-else-if="!gseaRes.rows.length" class="hint-line">{{ t('web.analysis.gsea_empty') }}</p>
          <template v-else>
            <div class="tbl-wrap">
              <table class="tbl">
                <thead><tr><th>term_name</th><th>NES</th><th>adj_p</th></tr></thead>
                <tbody>
                  <tr v-for="(r, i) in gseaRes.rows" :key="i">
                    <td>{{ r.Term ?? r.term_name ?? r.Name }}</td>
                    <td>{{ fmt(r.NES ?? r.nes, 3) }}</td>
                    <td>{{ fmt(r['FDR q-val'] ?? r.adj_p ?? r['FDR p-val'], 4) }}</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </template>
        </template>
      </details>
    </section>

    <!-- 模块 -->
    <section class="sec">
      <details class="fold" open>
        <summary class="sec-title fold-title">{{ t('web.analysis.sec_modules') }}</summary>
        <button class="bar-btn" :disabled="paramsStale || modsJob.running.value"
                @click="modsJob.run(api.rerunModules(ppi.directory_rel, ppi.prefix), { doneText: t('web.analysis.modules_done_text') })">
          {{ modsJob.running.value ? t('web.analysis.rerunning') : t('web.analysis.modules_rerun') }}
        </button>
        <p v-if="modsJob.running.value" class="hint-line">{{ modsJob.text.value || t('web.analysis.rerunning') }}</p>
        <p v-if="modsJob.error.value" class="notice err">{{ t('web.analysis.modules_done_text') }}:{{ modsJob.error.value }}</p>
        <p v-if="ppi.modules && !ppi.modules.summary.length" class="hint-line">
          {{ t('web.analysis.modules_empty') }}
        </p>
        <p v-if="ppi.modules?.summary.length" class="hint-line">
          {{ t('web.analysis.modules_summary_title') }}<template v-if="ppi.modules.hidden_single">
            · {{ t('web.analysis.modules_hidden', { n: ppi.modules.hidden_single }) }}</template>
        </p>
        <div v-if="ppi.modules?.summary.length" class="tbl-wrap">
          <table class="tbl">
            <thead><tr><th>module_id</th><th>size</th><th>avg_degree</th><th>genes_preview</th></tr></thead>
            <tbody>
              <tr v-for="r in ppi.modules.summary" :key="String(r.module_id)">
                <td>{{ r.module_id }}</td><td>{{ r.size }}</td><td>{{ r.avg_degree }}</td><td>{{ r.genes_preview }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <template v-if="ppi.modules?.enr_rows?.length">
          <h4 class="sub-sec">{{ t('web.analysis.mod_enr_title') }}</h4>
          <p class="hint-line">
            {{ t('web.analysis.mod_enr_stats', { total: ppi.modules.enr_total, mods: ppi.modules.enr_modules }) }}
            <a v-if="art('module_enrichment_csv')" class="bar-btn"
               :href="artifactsUrl(art('module_enrichment_csv') ?? '', true)">{{ t('web.analysis.dl_full_csv') }}</a>
          </p>
          <div class="tbl-wrap">
            <table class="tbl">
              <thead>
                <tr><th>module_id</th><th>source</th><th>term_id</th><th>term_name</th>
                    <th>adj_p</th><th>intersection</th><th>query_size</th></tr>
              </thead>
              <tbody>
                <tr v-for="(r, i) in ppi.modules.enr_rows" :key="i">
                  <td>{{ r.module_id }}</td><td>{{ r.source }}</td><td>{{ r.term_id }}</td>
                  <td>{{ r.term_name }}</td><td>{{ fmt(r.adjusted_p_value, 4) }}</td>
                  <td>{{ r.intersection_size }}</td><td>{{ r.query_size }}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </template>
        <p v-if="ppi.modules?.no_sig?.length" class="hint-line">
          {{ t('web.analysis.modules_no_sig', { names: ppi.modules.no_sig.join(', ') }) }}
        </p>
        <p v-if="ppi.modules?.failures?.length" class="stale">
          {{ t('web.analysis.modules_failures', { n: ppi.modules.failures.length }) }}
        </p>
      </details>
    </section>

    <!-- 下载区 -->
    <section class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">{{ t('web.analysis.dl_title', { n: downloadItems.length }) }}</summary>
        <button class="bar-btn" @click="zipWholeDirectory">{{ t('web.analysis.dl_whole', { name: `${dirName}.zip` }) }}</button>
        <ul class="dl-list">
          <li v-for="d in downloadItems" :key="d.rel" class="dl-row">
            <span class="dl-name">{{ d.rel }}</span>
            <a class="mini" :href="artifactsUrl(d.rel, true)">{{ t('web.common.download') }}</a>
          </li>
        </ul>
        <p class="hint-line">{{ t('web.analysis.dl_hint', { dir: ppi.directory_rel }) }}</p>
      </details>
    </section>
  </div>
</template>

<style scoped>
.mode-row {
  display: flex;
  align-items: baseline;
  gap: 16px;
  padding-bottom: 12px;
  border-bottom: 1px solid #e6e6e4;
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

.names-form {
  display: flex;
  align-items: flex-end;
  gap: 26px;
  margin-top: 26px;
}
.field {
  position: relative;
  width: clamp(180px, 17vw, 250px);
}
.field.genes {
  width: 100%;
  margin-top: 26px;
}
.field input,
.field textarea {
  width: 100%;
  box-sizing: border-box;
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid #e6e6e4;
  padding: 10px 2px;
  font-size: 15px;
  font-family: inherit;
  color: #1a1a1a;
  letter-spacing: 0.02em;
  border-radius: 0;
  resize: vertical;
}
.field input::placeholder,
.field textarea::placeholder {
  color: #cfcfcd;
  letter-spacing: 0.08em;
}
.field .line {
  position: absolute;
  left: 0;
  right: 0;
  bottom: 0;
  height: 1px;
  background: #1a1a1a;
  transform: scaleX(0);
  transition: transform 0.45s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.field:focus-within .line {
  transform: scaleX(1);
}
.cross {
  padding-bottom: 12px;
  color: #d7d7d5;
  font-size: 15px;
}

.opt-grid {
  display: flex;
  flex-wrap: wrap;
  gap: 22px 34px;
  margin-top: 26px;
}
.opt {
  display: flex;
  flex-direction: column;
  gap: 6px;
  min-width: 160px;
}
.opt-label {
  font-size: 10.5px;
  letter-spacing: 0.14em;
  color: #bcbcbc;
}
.opt input[type='number'],
.opt select {
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid #e6e6e4;
  padding: 6px 2px;
  font-family: inherit;
  font-size: 13.5px;
  color: #1a1a1a;
  border-radius: 0;
}
.opt input[type='range'] {
  accent-color: #1a1a1a;
  height: 2px;
  padding: 8px 0;
}

.build-row {
  margin-top: 34px;
  display: flex;
  align-items: center;
  gap: 16px;
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
.build:active:not(:disabled) {
  transform: scale(0.95);
}
.build:disabled {
  opacity: 0.4;
  cursor: default;
}
.build-hint {
  font-size: 11px;
  letter-spacing: 0.3em;
  color: #bcbcbc;
}

.prog {
  margin-top: 22px;
}
.stages {
  list-style: none;
  margin: 12px 0 0;
  padding: 0;
  display: flex;
  flex-wrap: wrap;
  gap: 8px 22px;
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
  margin-right: 7px;
  font-style: normal;
}

/* ---------------- 结果区 ---------------- */
.result {
  margin-top: 56px;
  padding-top: 26px;
  border-top: 1px solid #f0f0ee;
}
.sec {
  margin-top: 34px;
  padding-top: 20px;
  border-top: 1px solid #f0f0ee;
}
.sec-title {
  margin: 0 0 12px;
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0.14em;
  color: #1a1a1a;
}
.fold {
  border: none;
}
.fold-title {
  cursor: pointer;
  list-style: none;
}
.fold-title::before {
  content: '› ';
  color: #c2c2c2;
}
.fold[open] .fold-title::before {
  content: '⌄ ';
}
.tbl-wrap {
  max-height: 420px;
  overflow: auto;
  margin-top: 10px;
}
.tbl {
  width: 100%;
  border-collapse: collapse;
  font-size: 11.5px;
  color: #4a4a48;
}
.tbl th {
  text-align: left;
  font-weight: 500;
  letter-spacing: 0.08em;
  color: #a6a6a6;
  border-bottom: 1px solid #e6e6e4;
  padding: 6px 10px 6px 0;
  position: sticky;
  top: 0;
  background: #fdfdfd;
}
.tbl td {
  border-bottom: 1px solid #f4f4f2;
  padding: 6px 10px 6px 0;
  font-variant-numeric: tabular-nums;
}
.fig-frame {
  width: 100%;
  height: 640px;
  border: 1px solid #f0f0ee;
  margin-top: 12px;
  background: #ffffff;
}
.fig-img {
  max-width: 100%;
  margin-top: 12px;
  border: 1px solid #f0f0ee;
}
.bar-btn {
  border: none;
  background: none;
  font-family: inherit;
  font-size: 12px;
  letter-spacing: 0.08em;
  color: #1a1a1a;
  cursor: pointer;
  padding: 4px 0;
  margin-right: 18px;
  border-bottom: 1px solid #1a1a1a;
  transition: opacity 0.25s ease;
}
.bar-btn:disabled {
  cursor: default;
  border-bottom-color: #dedede;
  color: #c6c6c4;
}
.dl-list {
  list-style: none;
  margin: 10px 0 0;
  padding: 0;
}
.dl-row {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 14px;
  padding: 8px 2px;
  border-bottom: 1px solid #f4f4f2;
}
.dl-name {
  font-family: 'SF Mono', Consolas, monospace;
  font-size: 11px;
  color: #6e6e6c;
  word-break: break-all;
}
.notice.err {
  color: #b3261e;
}
</style>
