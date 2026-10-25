<script setup lang="ts">
/** PPI 分析(Data_Network 的 TS 迁移,阶段 3a)。
 * 产物直出:分析图表由后端落盘、经 /api/artifacts 伺服;交互区重跑(模块/GSEA)
 * 无状态化——从落盘 nodes/edges CSV 重建图。行为与 Streamlit 版一致。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import {
  api, artifactsUrl, type PpiAnalyzeBody, type PpiResult, type Row
} from '../api'
import { takePrefill } from '../store'
import { useJob } from '../composables/useJob'
import Stepper from '../components/Stepper.vue'

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

const STAGE_LABEL: Record<string, string> = {
  targets: '靶点获取与交集', string: 'STRING 互作查询', graph: '建图 · 拓扑 · Hub',
  enrichment: '功能富集', gsea: 'GSEA 排序富集', network: '邻近度 + RWR',
  modules: 'Louvain 模块检测', save: '结果落盘'
}
const STATE_ICON: Record<string, string> = { pending: '○', running: '◐', done: '●' }

const stepStates = computed<Record<string, string>>(() => {
  const map: Record<string, string> = {}
  for (const k of stages.value) map[k] = 'pending'
  const p = job.extra.value as { step?: string; state?: string }
  if (p?.step && p?.state && p.step in map) map[p.step] = p.state
  return map
})

function analyze(): void {
  if (mode.value === 'names' && (!drug.value.trim() || !disease.value.trim())) {
    flash('请填写药物成分与疾病表型(或切换到基因列表模式)')
    return
  }
  const body: PpiAnalyzeBody = {
    mode: mode.value, drug: drug.value, disease: disease.value, genes_text: genesText.value,
    score: score.value, network_type: networkType.value, rwr_expand: rwrExpand.value,
    hub_top_n: hubTopN.value, hub_method: hubMethod.value
  }
  void job.run(api.analyzePpi(body), { doneText: '报告已生成' })
}

/* ---------------- 结果 ---------------- */

const paramsStale = computed(() => {
  if (!ppi.value) return false
  const a = snapshot.value
  const b = ppi.value.params
  return a.length !== b.length || a.some((v, i) => String(v) !== String(b[i]))
})

const notice = ref('')
let noticeTimer = 0
function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 3200)
}

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
const VIEW_LABEL: Record<string, string> = {
  interactive: '交互式', static: '静态图', string: 'STRING 风格', cytoscape: 'Cytoscape 风格'
}

const artifacts = computed(() => ppi.value?.artifacts ?? {})
function art(key: string): string | null {
  return artifacts.value[key] ?? null
}

/* 邻近度 */
const prox = computed<Record<string, any> | null>(() =>
  (ppi.value?.proximity ?? null) as Record<string, any> | null)

/* RWR */
const rwr = computed<Record<string, Record<string, any>> | null>(() =>
  (ppi.value?.rwr ?? null) as Record<string, Record<string, any>> | null)
const rwrSide = ref<'disease_seeds' | 'drug_seeds'>('disease_seeds')
const rwrRes = computed<Record<string, any> | null>(() => rwr.value?.[rwrSide.value] ?? null)
const rwrTop = computed<Row[]>(() => {
  const res = rwrRes.value
  if (!res?.ok) return []
  const seeds: Row[] = ((res.seeds ?? []) as Row[]).map((r) => ({ ...r, is_seed: true }))
  const non: Row[] = ((res.top_nonseed ?? []) as Row[]).map((r) => ({ ...r, is_seed: false }))
  return [...seeds, ...non].sort((a, b) => Number(b.score) - Number(a.score))
})
const rwrNonSeed = computed(() => rwrTop.value.filter((r) => !r.is_seed))
const rwrSeeds = computed(() => rwrTop.value.filter((r) => r.is_seed))

/* 模块 / GSEA 重跑 */
const modsJob = useJob((result) => {
  if (ppi.value) {
    ppi.value = { ...ppi.value, modules: {
      summary: (result.summary ?? []) as Row[],
      failures: (result.failures ?? []) as unknown[],
      no_sig: (result.no_sig ?? []) as number[]
    } }
  }
  flash('模块检测已重跑完成')
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
  flash('GSEA 已重跑完成')
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
    'enrichment_bar_png', 'enrichment_bar_svg', 'enrichment_dot_png', 'enrichment_dot_svg']) {
    if (a[key]) out.push({ label: key, rel: a[key] })
  }
  return out
})

const dirName = computed(() =>
  (ppi.value?.directory_rel ?? '').split('/').pop() ?? '')

async function zipWholeDirectory(): Promise<void> {
  if (!dirName.value) return
  try {
    const blob = await api.zip('downloads', [dirName.value])
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${dirName.value}.zip`
    a.click()
    URL.revokeObjectURL(url)
  } catch {
    flash('整包下载失败')
  }
}

const fmt = (v: unknown, digits = 3): string =>
  typeof v === 'number' ? v.toFixed(digits) : (v == null ? '—' : String(v))

onMounted(() => {
  // 主页四库检索完成后的预填跳转
  const pf = takePrefill()
  if (pf) {
    mode.value = 'names'
    drug.value = pf.drug
    disease.value = pf.disease
  }
})

onBeforeUnmount(() => window.clearTimeout(noticeTimer))
</script>

<template>
  <p class="kicker s1">PPI&nbsp;ANALYSIS</p>
  <h1 class="title s2">网络药理学分析</h1>
  <p class="sub s3">成分&nbsp;—&nbsp;靶点&nbsp;—&nbsp;网络&nbsp;—&nbsp;富集&nbsp;—&nbsp;模块</p>

  <!-- ================= 表单 ================= -->
  <div class="panel wide s4">
    <div class="mode-row">
      <span class="ctl-label">输入模式</span>
      <span class="lv-group">
        <button class="lv" :class="{ on: mode === 'names' }" @click="mode = 'names'">药物 + 疾病名称</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: mode === 'genes' }" @click="mode = 'genes'">基因列表</button>
      </span>
    </div>

    <div v-if="mode === 'names'" class="names-form">
      <label class="field">
        <input v-model="drug" spellcheck="false" placeholder="药物成分" />
        <i class="line"></i>
      </label>
      <span class="cross">×</span>
      <label class="field">
        <input v-model="disease" spellcheck="false" placeholder="疾病表型" />
        <i class="line"></i>
      </label>
    </div>
    <label v-else class="field genes">
      <textarea v-model="genesText" rows="2" spellcheck="false"
                placeholder="基因符号,逗号/分号/空格分隔" />
      <i class="line"></i>
    </label>

    <div class="opt-grid">
      <label class="opt">
        <span class="opt-label">STRING 置信度 {{ score.toFixed(2) }}</span>
        <input v-model.number="score" type="range" min="0.15" max="0.9" step="0.05" />
      </label>
      <label class="opt">
        <span class="opt-label">网络类型</span>
        <select v-model="networkType">
          <option value="full">功能关联网络 (full)</option>
          <option value="physical">物理结合网络 (physical)</option>
        </select>
      </label>
      <label class="opt">
        <span class="opt-label">RWR 扩展靶点数(0 = 关闭)</span>
        <Stepper v-model="rwrExpand" :min="0" :max="100" />
      </label>
      <label class="opt">
        <span class="opt-label">Hub 数量(5–20)</span>
        <Stepper v-model="hubTopN" :min="5" :max="20" />
      </label>
      <label class="opt">
        <span class="opt-label">Hub 排名依据</span>
        <select v-model="hubMethod">
          <option value="degree">度数(文献惯例)</option>
          <option value="mcc">MCC(CytoHubba)</option>
        </select>
      </label>
    </div>

    <div class="build-row">
      <button class="build" :disabled="job.running.value || modsJob.running.value || gseaJob.running.value"
              title="生成报告" aria-label="生成报告" @click="analyze">
        <svg viewBox="0 0 24 14" width="26" height="15">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <span class="build-hint">生成报告</span>
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
          <i class="stage-icon">{{ STATE_ICON[stepStates[k]] }}</i>{{ STAGE_LABEL[k] }}
        </li>
      </ul>
    </div>
    <p v-if="job.error.value" class="notice err">分析失败:{{ job.error.value }}(完整堆栈见日志页)</p>
  </div>

  <!-- ================= 结果 ================= -->
  <div v-if="ppi" class="panel wide result">
    <p v-if="paramsStale" class="stale">
      表单参数已变化 —— 下方结果属于上一次参数;重新点击"生成报告"刷新。
    </p>

    <div class="metrics s5">
      <div class="metric"><b>{{ ppi.genes_count }}</b><span>输入靶点</span></div>
      <div class="metric"><b>{{ ppi.nodes }}</b><span>网络节点</span></div>
      <div class="metric"><b>{{ ppi.edges }}</b><span>边</span></div>
      <div class="metric"><b>{{ ppi.hubs }}</b><span>Hub 基因</span></div>
    </div>

    <p v-if="ppi.expansion" class="hint-line">
      RWR 扩展:种子 {{ ppi.expansion.seeds }} + 扩展 {{ ppi.expansion.expanded }} =
      {{ ppi.expansion.seeds + ppi.expansion.expanded }} 个查询基因(α=0.85,前缀 {{ ppi.prefix }})
    </p>
    <p v-if="!ppi.genes_count" class="stale">无共同靶点 —— 请检查药物/疾病名称。</p>
    <p v-if="!ppi.edges" class="stale">
      STRING 未返回高于置信度阈值的互作{{ ppi.views ? '' : '' }}
      (可尝试降低置信度阈值或改用功能关联网络)。
    </p>

    <!-- Hub 表 -->
    <section v-if="ppi.hub_table.length" class="sec">
      <h3 class="sec-title">
        {{ ppi.hub_rule.startsWith('top')
          ? `Hub 基因 Top ${ppi.hub_table.length}(${ppi.hub_method_eff} 排名)` : 'Hub 基因' }}
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
      <h3 class="sec-title">网络图</h3>
      <div class="lv-group">
        <button v-for="v in ppi.views" :key="v" class="lv"
                :class="{ on: effectiveView === v }" @click="setView(v)">{{ VIEW_LABEL[v] }}</button>
      </div>
      <p v-if="effectiveView === 'string'" class="hint-line">STRING 数据库原生风格:边粗细 = 证据强度,节点颜色 = 度数。</p>
      <p v-if="effectiveView === 'cytoscape'" class="hint-line">Cytoscape 风格:度数降序同心圆排布,颜色映射模块归属。</p>
      <iframe v-if="effectiveView === 'interactive'"
              class="fig-frame" :src="artifactsUrl(art('html') ?? '')" />
      <img v-else-if="effectiveView === 'static'"
           class="fig-img" :src="artifactsUrl(art('png') ?? '')" alt="static network" />
      <img v-else-if="effectiveView === 'string'"
           class="fig-img" :src="artifactsUrl(art('string_png') ?? '')" alt="string style" />
      <img v-else class="fig-img" :src="artifactsUrl(art('cyto_png') ?? '')" alt="cytoscape style" />
      <p class="hint-line">四种视图均已落盘,可在下方下载区取用。</p>
    </section>

    <!-- 核心靶点 -->
    <section v-if="art('core_png')" class="sec">
      <h3 class="sec-title">核心靶点子网络(MCC Top 10)</h3>
      <img class="fig-img" :src="artifactsUrl(art('core_png') ?? '')" alt="core targets" />
    </section>

    <!-- 邻近度 -->
    <section v-if="prox && ppi.params[0] === 'names'" class="sec">
      <h3 class="sec-title">网络邻近度检验</h3>
      <p class="hint-line">
        度匹配零分布 {{ prox.n_perm }} 次置换 · 种子 {{ prox.seed ?? '-' }} · 度分箱
        {{ 100 }}
      </p>
      <template v-if="prox.error === 'not_built'">
        <p class="stale">背景互作组未构建 —— 请先到数据管理页构建。</p>
      </template>
      <template v-else-if="['low_coverage_drug', 'low_coverage_disease', 'empty_after_filter'].includes(prox.error)">
        <p class="stale">
          靶点落入背景网络的比例过低(药物 {{ prox.coverage?.drug_in_lcc ?? 0 }}/{{ prox.coverage?.drug_total ?? 0 }},
          疾病 {{ prox.coverage?.disease_in_lcc ?? 0 }}/{{ prox.coverage?.disease_total ?? 0 }}),检验不可靠。
        </p>
      </template>
      <template v-else-if="prox.error">
        <p class="stale">邻近度检验失败:{{ prox.error }}</p>
      </template>
      <template v-else-if="prox.ok">
        <div class="metrics three">
          <div class="metric"><b>{{ fmt(prox.d_c) }}</b><span>d_C(平均最短距离)</span></div>
          <div class="metric"><b>{{ fmt(prox.z) }}</b><span>Z-score</span></div>
          <div class="metric"><b>{{ fmt(prox.p_empirical, 4) }}</b><span>p(经验)</span></div>
        </div>
        <p class="hint-line" :class="{ good: typeof prox.z === 'number' && prox.z < 0 }">
          {{ typeof prox.z === 'number'
            ? (prox.z < 0 ? '药物靶点与疾病蛋白的网络距离显著近于随机期望 —— 支持二者存在作用关系。'
               : '网络距离未优于随机期望 —— 不支持存在特异性作用关系。')
            : '零分布退化,Z 无意义(两靶点集过近)。' }}
        </p>
        <p class="hint-line">
          背景网络 LCC {{ prox.lcc_nodes }} 节点 ·
          覆盖:药物 {{ prox.coverage?.drug_in_lcc ?? 0 }}/{{ prox.coverage?.drug_total ?? 0 }},
          疾病 {{ prox.coverage?.disease_in_lcc ?? 0 }}/{{ prox.coverage?.disease_total ?? 0 }}
          {{ prox.cached ? '· 缓存命中' : '' }}
        </p>
        <template v-if="prox.d_c_excl != null">
          <p class="hint-line">去交集版本(剔除 {{ prox.n_excl ?? 0 }} 个共同靶点,排除"纯属直接作用"的解释):</p>
          <div class="metrics three">
            <div class="metric"><b>{{ fmt(prox.d_c_excl) }}</b><span>d_C(excl)</span></div>
            <div class="metric"><b>{{ fmt(prox.z_excl) }}</b><span>Z(excl)</span></div>
            <div class="metric"><b>{{ fmt(prox.p_excl, 4) }}</b><span>p(excl)</span></div>
          </div>
        </template>
      </template>
    </section>

    <!-- RWR -->
    <section v-if="rwr && ppi.params[0] === 'names'" class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">靶点优先级排序(RWR)</summary>
        <p class="hint-line">α = 0.85 重启概率;单侧切换只查看已算结果,不重算。</p>
        <div class="lv-group">
          <button class="lv" :class="{ on: rwrSide === 'disease_seeds' }"
                  @click="rwrSide = 'disease_seeds'">疾病种子</button>
          <i class="lv-sep">/</i>
          <button class="lv" :class="{ on: rwrSide === 'drug_seeds' }"
                  @click="rwrSide = 'drug_seeds'">药物种子</button>
        </div>
        <template v-if="rwrRes">
          <p v-if="rwrRes.error === 'not_built'" class="stale">背景互作组未构建。</p>
          <p v-else-if="rwrRes.error" class="stale">RWR 失败:{{ rwrRes.error }}</p>
          <template v-else-if="rwrRes.ok">
            <p class="hint-line">
              迭代 {{ rwrRes.n_iter }} 次,{{ rwrRes.converged ? '已收敛' : '未收敛' }} ·
              种子入网 {{ rwrRes.coverage?.seed_in_lcc ?? 0 }}/{{ rwrRes.coverage?.seed_total ?? 0 }}
              {{ rwrRes.cached ? '· 缓存命中' : '' }}
            </p>
            <p v-if="rwrRes.enrichment?.z != null" class="hint-line"
               :class="{ good: rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05 }">
              度匹配校正:标注集 {{ rwrRes.enrichment.n_annotated }} 个,
              Z = {{ fmt(rwrRes.enrichment.z) }},
              p = {{ fmt(rwrRes.enrichment.p, 4) }}
              {{ rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05 ? '(显著富集)' : '(不显著)' }}
            </p>
            <div v-if="rwrNonSeed.length" class="tbl-wrap">
              <table class="tbl">
                <thead><tr><th>gene</th><th>score</th><th v-if="'annotated' in rwrNonSeed[0]">标注集</th></tr></thead>
                <tbody>
                  <tr v-for="r in rwrNonSeed" :key="String(r.gene)">
                    <td>{{ r.gene }}</td><td>{{ fmt(r.score, 6) }}</td>
                    <td v-if="'annotated' in r">{{ r.annotated ? '✓' : '' }}</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p v-else class="hint-line">无排序结果。</p>
            <details v-if="rwrSeeds.length" class="fold">
              <summary class="hint-line">种子基因({{ rwrSeeds.length }})</summary>
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
      <h3 class="sec-title">功能富集</h3>
      <p v-if="ppi.enrichment.reason === 'both_sources_down' || ppi.enrichment.reason === 'gprofiler_down_no_fallback'"
         class="stale">
        富集双源均不可用,结果不可信 —— 稍后自动重试或手动重跑主分析。
      </p>
      <p v-else-if="!ppi.enrichment.rows.length" class="hint-line">无显著富集项。</p>
      <template v-else>
        <p class="hint-line">
          数据源 {{ ppi.enrichment.source }} · {{ ppi.enrichment.rows.length }} 项显著通路
        </p>
        <div class="tbl-wrap">
          <table class="tbl">
            <thead><tr><th>来源</th><th>term_id</th><th>term_name</th><th>adj_p</th><th>交集</th><th>查询集</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in ppi.enrichment.rows" :key="i">
                <td>{{ r.source }}</td><td>{{ r.term_id }}</td><td>{{ r.term_name }}</td>
                <td>{{ fmt(r.adjusted_p_value, 4) }}</td><td>{{ r.intersection_size }}</td><td>{{ r.query_size }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <template v-if="art('enrichment_bar_png')">
          <p class="hint-line">显著通路柱状图(Top 15):</p>
          <img class="fig-img" :src="artifactsUrl(art('enrichment_bar_png') ?? '')" alt="enrichment bar" />
        </template>
        <template v-if="art('enrichment_dot_png')">
          <p class="hint-line">气泡图:</p>
          <img class="fig-img" :src="artifactsUrl(art('enrichment_dot_png') ?? '')" alt="enrichment dot" />
        </template>
        <template v-if="art('dtp_png')">
          <p class="hint-line">药物-靶点-通路网络图(KEGG):</p>
          <img class="fig-img" :src="artifactsUrl(art('dtp_png') ?? '')" alt="drug target pathway" />
        </template>
      </template>
    </section>

    <!-- GSEA -->
    <section v-if="ppi.params[0] === 'names'" class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">GSEA 排序富集</summary>
        <div class="lv-group">
          <button class="lv" :class="{ on: gseaSide === 'drug' }" @click="gseaSide = 'drug'">药物侧</button>
          <i class="lv-sep">/</i>
          <button class="lv" :class="{ on: gseaSide === 'disease' }" @click="gseaSide = 'disease'">疾病侧</button>
        </div>
        <button class="bar-btn" :disabled="gseaDone || paramsStale || gseaJob.running.value"
                @click="gseaSide && gseaJob.run(api.rerunGsea(gseaSide, drug, disease, ppi.prefix), { doneText: 'GSEA 完成' })">
          {{ gseaJob.running.value ? '重跑中…' : (gseaDone ? '该侧已有结果' : '重跑该侧') }}
        </button>
        <p v-if="gseaJob.running.value" class="hint-line">{{ gseaJob.text.value || '计算中…' }}</p>
        <p v-if="gseaJob.error.value" class="notice err">GSEA 失败:{{ gseaJob.error.value }}</p>
        <template v-if="gseaRes">
          <p v-if="gseaRes.reason === 'insufficient_input'" class="stale">可用基因数不足,跳过。</p>
          <p v-else-if="gseaRes.reason === 'request_failed'" class="stale">请求失败 —— 可点击"重跑该侧"重试。</p>
          <p v-else-if="!gseaRes.rows.length" class="hint-line">无通过 FDR 筛选的通路。</p>
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
        <summary class="sec-title fold-title">Louvain 模块检测</summary>
        <button class="bar-btn" :disabled="paramsStale || modsJob.running.value"
                @click="modsJob.run(api.rerunModules(ppi.directory_rel, ppi.prefix), { doneText: '模块重跑完成' })">
          {{ modsJob.running.value ? '重跑中…' : '重新运行模块检测与富集' }}
        </button>
        <p v-if="modsJob.running.value" class="hint-line">{{ modsJob.text.value || '计算中…' }}</p>
        <p v-if="modsJob.error.value" class="notice err">模块重跑失败:{{ modsJob.error.value }}</p>
        <p v-if="ppi.modules && !ppi.modules.summary.length" class="hint-line">
          无规模达标的功能模块(网络太小或全孤立)。
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
        <p v-if="ppi.modules?.no_sig?.length" class="hint-line">
          无显著通路的模块:{{ ppi.modules.no_sig.join(', ') }}
        </p>
        <p v-if="ppi.modules?.failures?.length" class="stale">
          富集失败模块:{{ ppi.modules.failures.length }} 个(详情见日志)。
        </p>
      </details>
    </section>

    <!-- 下载区 -->
    <section class="sec">
      <details class="fold">
        <summary class="sec-title fold-title">下载产物(共 {{ downloadItems.length }} 项 + 整包)</summary>
        <button class="bar-btn" @click="zipWholeDirectory">整包下载 ({{ dirName }}.zip)</button>
        <ul class="dl-list">
          <li v-for="d in downloadItems" :key="d.rel" class="dl-row">
            <span class="dl-name">{{ d.rel }}</span>
            <a class="mini" :href="artifactsUrl(d.rel, true)">下载</a>
          </li>
        </ul>
        <p class="hint-line">全部产物位于 {{ ppi.directory_rel }},也可在数据管理页整目录下载。</p>
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
.stale {
  margin: 10px 0 0;
  font-size: 12px;
  letter-spacing: 0.05em;
  color: #b3261e;
}
.metrics {
  display: flex;
  gap: 56px;
  justify-content: center;
  margin: 8px 0 6px;
}
.metrics.three {
  gap: 40px;
  justify-content: flex-start;
}
.metric {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 2px;
}
.metric b {
  font-size: 30px;
  font-weight: 300;
  font-variant-numeric: tabular-nums;
  color: #1a1a1a;
}
.metric span {
  font-size: 10.5px;
  letter-spacing: 0.2em;
  color: #bcbcbc;
}
.hint-line {
  margin: 10px 0 0;
  font-size: 11.5px;
  letter-spacing: 0.05em;
  color: #8a8a88;
}
.hint-line.good {
  color: #1a1a1a;
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
