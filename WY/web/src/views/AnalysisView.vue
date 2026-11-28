<script setup lang="ts">
/** PPI 分析(Data_Network 的 TS 迁移,阶段 3a)。
 * 产物直出:分析图表由后端落盘、经 /api/artifacts 伺服;交互区重跑(模块/GSEA)
 * 无状态化——从落盘 nodes/edges CSV 重建图。行为与 Streamlit 版一致。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import {
  api, artifactsUrl, type PpiResult, type Row
} from '../api'
import { takePrefill } from '../store'
import { useNotice } from '../composables/useNotice'
import { useAnalysisJob } from '../composables/useAnalysisJob'
import Stepper from '../components/Stepper.vue'
import AnalysisDownloadSection from '../components/AnalysisDownloadSection.vue'
import AnalysisGseaSection from '../components/AnalysisGseaSection.vue'
import AnalysisModulesSection from '../components/AnalysisModulesSection.vue'
import AnalysisNetworkViewer from '../components/AnalysisNetworkViewer.vue'
import AnalysisProximitySection from '../components/AnalysisProximitySection.vue'
import AnalysisReachCmpSection from '../components/AnalysisReachCmpSection.vue'
import AnalysisRwrSection from '../components/AnalysisRwrSection.vue'
import ConsensusPanel from '../components/ConsensusPanel.vue'
import { fmt } from '../utils/format'

const { t } = useI18n()

/* ---------------- 表单 ---------------- */

/* 已缓存查询名下拉(与协同页同源同款):本地 merged 缓存的药/病查询名,
 * 只读本地不联网;拉取失败输入框照常手输(SynergyView 同款降级) */
const drugOptions = ref<string[]>([])
const diseaseOptions = ref<string[]>([])

const mode = ref<'names' | 'genes'>('names')
const drug = ref('')
const disease = ref('')
const genesText = ref('')
const score = ref(0.45)
const networkType = ref('full')
const rwrExpand = ref(0)
const hubTopN = ref(10)
const hubMethod = ref('degree')

const snapshot = computed<(string | number)[]>(() => [
  mode.value,
  mode.value === 'names' ? drug.value.trim() : '',
  mode.value === 'names' ? disease.value.trim() : '',
  genesText.value.trim(),
  score.value, networkType.value, rwrExpand.value, hubTopN.value, hubMethod.value
])

const { notice, flash } = useNotice(3200)

/* 主作业(ppi 数据源/job/七阶段显示/hub 表可达性/paramsStale/analyze/resume)
 * 见 useAnalysisJob;表单状态留壳层——一致性面板 fill 联动写表单字段,datalist
 * 与模板绑定都在此处,拆走接线大于收益 */
const {
  ppi, paramsStale, job, reachFilter, hasReach, reachTier, reachWord, reachTitle,
  hubCols, hubRowsShown, stages, stageLabels, STATE_ICON, stepStates,
  analyze, resumeAnalyzeJob
} = useAnalysisJob({
  mode, drug, disease, genesText, score, networkType, rwrExpand, hubTopN, hubMethod,
  snapshot, flash
})

const artifacts = computed(() => ppi.value?.artifacts ?? {})
function art(key: string): string | null {
  return artifacts.value[key] ?? null
}

/* 一致性面板强制刷新信号:主分析 job 的 jobId 每次提交唯一——同名药病对重跑
 * 时运行目录不变,directory_rel 观察不到;GSEA/模块重跑走独立端点不碰主 jobId */
const runKey = computed(() => job.jobId.value)

/* 邻近度(PpiResult.proximity 已是 ProximityResult,无需再转 any) */
const prox = computed(() => ppi.value?.proximity ?? null)

/* 暴露约束前后排序对比摘要(W7):gate=false → 结论行;true → artifacts 出图 */
const reachCmp = computed(() => ppi.value?.reach_compare ?? null)

/* RWR(侧切换与行拆分在 AnalysisRwrSection) */
const rwr = computed(() => ppi.value?.rwr ?? null)

/* 模块 / GSEA 重跑:完成事件由组件 emit,此处合并回 PpiResult(单一结果对象)
 * 并 flash;api.rerunGsea/rerunModules 调用已移入组件 script */
function onGseaDone(result: Record<string, unknown>, side: 'drug' | 'disease'): void {
  if (!ppi.value) return
  const g = { ...ppi.value.gsea }
  g[side] = {
    reason: String(result.reason ?? ''),
    stats: (result.stats ?? {}) as Row,
    rows: (result.rows ?? []) as Row[]
  }
  ppi.value = { ...ppi.value, gsea: g }
  flash(t('web.analysis.gsea_done_flash'))
}
function onModulesDone(modules: NonNullable<PpiResult['modules']>): void {
  if (!ppi.value) return
  ppi.value = { ...ppi.value, modules }
  flash(t('web.analysis.modules_done_flash'))
}
/* 重跑忙闲回传(组件 emit):「开始分析」按钮防重入——重跑中禁点 */
const gseaRerunning = ref(false)
const modsRerunning = ref(false)

/* 多方法一致性面板(模块C)联动:
 * - saved:一致性产物是在分析完成后另起的任务写的,本次 ppi.files 快照里没有
 *   ——按后端 files_written 补进下载列表(目录 rel 与 _serialize_ppi 同构);
 * - fill(下游打通 v1):一致子集填进手动模式基因框,参数过期标记自然亮起,
 *   用户自己点「开始分析」重跑——严守方案 §4.4「默认不替换」红线。 */
function consensusSaved(names: string[]): void {
  const dir = ppi.value?.directory_rel
  if (!dir || !ppi.value) return
  const known = new Set((ppi.value.files ?? []).map((f) => f.name))
  const add = names.filter((n) => n && !known.has(n))
    .map((n) => ({ name: n, rel: `${dir}/${n}` }))
  if (add.length) ppi.value = { ...ppi.value, files: [...(ppi.value.files ?? []), ...add] }
}

function fillFromConsensus(genes: string[]): void {
  if (!genes.length) return
  mode.value = 'genes'
  genesText.value = genes.join(', ')
  flash(t('web.analysis.consensus_filled', { count: genes.length }))
}

onMounted(() => {
  // 主页四库检索完成后的预填跳转
  const pf = takePrefill()
  if (pf) {
    mode.value = 'names'
    drug.value = pf.drug
    disease.value = pf.disease
  }
  void resumeAnalyzeJob()
  void (async () => {
    try {
      const lib = await api.synergyLibrary()
      drugOptions.value = lib.drugs
      diseaseOptions.value = lib.diseases
    } catch {
      /* 本地库清单不可达时输入框照常手输 */
    }
  })()
})
</script>

<template>
  <p class="kicker s1">{{ t('web.analysis.kicker') }}</p>
  <h1 class="title s2">{{ t('web.analysis.title') }}</h1>
  <p class="sub s3">{{ t('web.analysis.sub') }}</p>

  <!-- ================= 表单 ================= -->
  <div class="panel wide s4">
    <div class="xlinks">
      <a class="bar-btn" href="#/artifacts">{{ t('web.nav.artifacts_link') }}</a>
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
      <!-- 已缓存查询名下拉(2026-11-28 负责人定稿):与协同页同款 datalist + 自由手输,
           数据源 /api/synergy/library(list_known_targets,merged 缓存命中名)——
           选中即缓存命中,分析免重新联网,一致性面板种子解析也秒回 -->
      <datalist id="analysis-drug-options">
        <option v-for="d in drugOptions" :key="d" :value="d"></option>
      </datalist>
      <datalist id="analysis-disease-options">
        <option v-for="d in diseaseOptions" :key="d" :value="d"></option>
      </datalist>
      <label class="field">
        <input v-model="drug" spellcheck="false" list="analysis-drug-options" :placeholder="t('web.common.ph_drug')" />
        <i class="line"></i>
      </label>
      <span class="cross">×</span>
      <label class="field">
        <input v-model="disease" spellcheck="false" list="analysis-disease-options" :placeholder="t('web.common.ph_disease')" />
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
      <button class="build" :disabled="job.running.value || gseaRerunning || modsRerunning"
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
    <!-- 衔接说明(2026-11-28 实测 Ferulic acid×T2DM 后定稿):0 交集时下方邻近度/
         RWR/一致性仍有结果,不解释会让用户误以为"没网还能算"是 bug——它们走的是
         背景互作组 × 双方全量靶点的另一条数据线,与交集网络无关 -->
    <p v-if="!ppi.genes_count" class="hint-line">{{ t('web.analysis.empty_overlap_hint') }}</p>

    <!-- Hub 表(可达性标注列 W7:names 模式懒取;过滤默认"全部"保一期) -->
    <section v-if="ppi.hub_table.length" class="sec">
      <h3 class="sec-title">
        {{ ppi.hub_rule.startsWith('top')
          ? t('web.analysis.hub_top_title', { n: ppi.hub_table.length, method: ppi.hub_method_eff })
          : t('web.analysis.hub_title') }}
        <select v-if="hasReach" v-model="reachFilter" class="reach-filter"
                :title="t('web.reach.filter')" style="margin-left:12px">
          <option value="all">{{ t('web.reach.filter_all') }}</option>
          <option value="partial_up">{{ t('web.reach.filter_partial_up') }}</option>
          <option value="reachable">{{ t('web.reach.filter_reachable') }}</option>
          <option value="judged">{{ t('web.reach.filter_judged') }}</option>
        </select>
      </h3>
      <div class="tbl-wrap">
        <table class="tbl">
          <thead>
            <tr>
              <th v-for="k in hubCols" :key="k">{{ k === '__reach' ? t('web.reach.col') : k }}</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(r, i) in hubRowsShown" :key="i">
              <td v-for="(v, k) in r" :key="k">{{ v }}</td>
              <td v-if="hasReach">
                <span class="reach-cell" :class="`reach-${reachTier(r.gene)}`"
                      :title="reachTitle(r.gene)">
                  {{ reachWord(r.gene) }}
                </span>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </section>

    <!-- 网络图:视图切换 + 核心靶点 -->
    <AnalysisNetworkViewer :views="ppi.views" :artifacts="artifacts" />

    <!-- 邻近度(有结果且 names 模式才渲染;空网络时 proximity=null 整块不出现) -->
    <template v-if="prox && ppi.params[0] === 'names'">
      <AnalysisProximitySection :prox="prox" />
    </template>

    <!-- 暴露约束前后排序对比(W7 招牌图) -->
    <AnalysisReachCmpSection v-if="reachCmp" :reach-cmp="reachCmp" :artifacts="artifacts" />

    <!-- RWR(有结果且 names 模式才渲染;空网络时 rwr={} 整块不出现) -->
    <template v-if="rwr && Object.keys(rwr).length && ppi.params[0] === 'names'">
      <AnalysisRwrSection :rwr="rwr" />
    </template>

    <!-- 多方法一致性(模块C):独立组件,种子来源随本次分析(names=疾病靶点全集/genes=基因列表),
         产物落运行目录自动进下载列表;填入按钮走 fill 事件写表单,默认不替换 -->
    <ConsensusPanel v-if="ppi" :ppi="ppi" :mode="mode" :drug="drug" :disease="disease"
                    :genes-text="genesText" :run-key="runKey"
                    @fill="fillFromConsensus" @saved="consensusSaved" />

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

    <!-- GSEA(names 模式渲染条件在壳层) -->
    <template v-if="ppi.params[0] === 'names'">
      <AnalysisGseaSection :gsea="ppi.gsea" :artifacts="artifacts" :params-stale="paramsStale"
                           :drug="drug" :disease="disease" :prefix="ppi.prefix"
                           @done="onGseaDone" @running="v => gseaRerunning = v" />
    </template>

    <!-- 模块 -->
    <AnalysisModulesSection :modules="ppi.modules" :params-stale="paramsStale"
                            :directory-rel="ppi.directory_rel" :prefix="ppi.prefix"
                            @done="onModulesDone" @running="v => modsRerunning = v" />

    <!-- 下载区 -->
    <AnalysisDownloadSection :files="ppi.files" :directory-rel="ppi.directory_rel" :flash="flash" />
  </div>
</template>

<style scoped>
.mode-row {
  display: flex;
  align-items: baseline;
  gap: 16px;
  padding-bottom: 12px;
  border-bottom: 1px solid var(--c-line);
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
  border-bottom: 1px solid var(--c-line);
  padding: 10px 2px;
  font-size: 15px;
  font-family: inherit;
  color: var(--c-ink);
  letter-spacing: 0.02em;
  border-radius: 0;
  resize: vertical;
}
.field input::placeholder,
.field textarea::placeholder {
  color: var(--c-faint);
  letter-spacing: 0.08em;
}
.field .line {
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
  color: var(--c-ink3);
}
.opt input[type='number'],
.opt select {
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid var(--c-line);
  padding: 6px 2px;
  font-family: inherit;
  font-size: 13.5px;
  color: var(--c-ink);
  border-radius: 0;
}
.opt input[type='range'] {
  accent-color: var(--c-ink);
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
  background: var(--c-accent);
  color: var(--c-on-accent);
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
  color: var(--c-ink3);
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
  color: var(--c-ink4);
}
.stage.running {
  color: var(--c-ink);
}
.stage.done {
  color: var(--c-ink2);
}
.stage-icon {
  margin-right: 7px;
  font-style: normal;
}

/* ---------------- 结果区 ---------------- */
.result {
  margin-top: 56px;
  padding-top: 26px;
  border-top: 1px solid var(--c-line-soft);
}
.sec {
  margin-top: 34px;
  padding-top: 20px;
  border-top: 1px solid var(--c-line-soft);
}
.sec-title {
  margin: 0 0 12px;
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0.14em;
  color: var(--c-ink);
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
  color: var(--c-body);
}
.tbl th {
  text-align: left;
  font-weight: 500;
  letter-spacing: 0.08em;
  color: var(--c-sub);
  border-bottom: 1px solid var(--c-line);
  padding: 6px 10px 6px 0;
  position: sticky;
  top: 0;
  background: var(--c-paper);
}
.tbl td {
  border-bottom: 1px solid var(--c-rowline);
  padding: 6px 10px 6px 0;
  font-variant-numeric: tabular-nums;
}
.fig-frame {
  width: 100%;
  height: 640px;
  border: 1px solid var(--c-line-soft);
  margin-top: 12px;
  background: var(--c-panel);
}
.bar-btn {
  border: none;
  background: none;
  font-family: inherit;
  font-size: 12px;
  letter-spacing: 0.08em;
  color: var(--c-ink);
  cursor: pointer;
  padding: 4px 0;
  margin-right: 18px;
  border-bottom: 1px solid var(--c-accent);
  transition: opacity 0.25s ease;
}
.bar-btn:disabled {
  cursor: default;
  border-bottom-color: var(--c-ghost);
  color: var(--c-ink4);
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
  border-bottom: 1px solid var(--c-rowline);
}
/* 无友好名映射的文件:文件名暗淡斜体显示(可见但明确"未登记",提示维护者补映射) */
.dl-name.unknown {
  color: var(--c-ink3);
  font-style: italic;
}
/* dl-list/dl-main/dl-desc/dl-name 已提升为全局样式(styles/screen.css) */
.notice.err {
  color: var(--c-danger);
}
</style>
