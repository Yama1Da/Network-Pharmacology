<script setup lang="ts">
/** 疾病模块多方法一致性面板(二期方案 §4.4,模块C):RWR/DIAMOnD/近邻扩张
 * 三方法同一输入、同一规模 N 独立扩张,Borda 平均秩共识,被 ≥k 法捞入标记一致。
 * 静态 import(不懒加载):本地服务下懒加载只引入首开请求瀑布,抽组件的目的
 * 是 AnalysisView 瘦身(1056 行)而非省 bundle。
 * 热力配色 viridis(科研热图事实标准:感知均匀、色盲安全、黑白打印亮度单调),
 * 秩 1 = 最深端;「–」= 未被该方法捞入(非成员,不染色不误导)。色值是数据
 * 编码不属于界面皮肤,深浅主题下同一套(SynergyResultPanel 的 var(--c-*) 纪律
 * 只约束 chrome)。
 * 下游打通 v1:仅「将一致子集填入分析输入」按钮(emit fill,父组件写表单、
 * 用户自己重跑),严守方案 §4.4「默认不替换」红线;消融对比归 benchmark/CLI
 * 层,前端不做。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type ModulesConsensusResult } from '../api'
import { useJob } from '../composables/useJob'
import { fmt } from '../utils/format'

const props = defineProps<{
  ppi: { directory_rel: string; prefix: string }
  mode: 'names' | 'genes'
  drug: string
  disease: string
  genesText: string
}>()
const emit = defineEmits<{
  (e: 'fill', genes: string[]): void
  (e: 'saved', files: string[]): void
}>()
const { t } = useI18n()

const nRef = ref<number | ''>(100)
const kRef = ref<number | ''>(2)
const permsRef = ref<number | ''>(0)
const result = ref<ModulesConsensusResult | null>(null)

/* 每次新分析(= 新运行目录)都清掉上次一致性结果,强制用户重跑:一致性产物
 * 落在运行目录、经 ppi.files 快照进下载区——残留旧图会让用户误以为已有结果
 * 而跳过重跑,新目录便缺一致性产物(2026-11-28 用户定稿)。重跑模块/GSEA
 * 不改运行目录,不触发本钩子。 */
watch(() => props.ppi?.directory_rel, () => {
  result.value = null
})

/* 参数校验(2026-11-28 负责人质询后加固):越界值不许静默钳位执行——实时红标
 * + 明示范围 + 拒绝提交;空框在提交前回写默认值,让"按默认跑"可见而非静默。
 * 后端 400 范围校验保留为最后一道闸,正常路径到不了那里。 */
const effN = computed(() => (nRef.value === '' ? 100 : Math.round(Number(nRef.value))))
const effK = computed(() => (kRef.value === '' ? 2 : Math.round(Number(kRef.value))))
const effPerms = computed(() => (permsRef.value === '' ? 0 : Math.round(Number(permsRef.value))))
const invalidN = computed(() => nRef.value !== '' && !(effN.value >= 10 && effN.value <= 300))
const invalidK = computed(() => kRef.value !== '' && !(effK.value >= 1 && effK.value <= 3))
const invalidPerms = computed(() => permsRef.value !== '' && !(effPerms.value >= 0 && effPerms.value <= 200))
const anyInvalid = computed(() => invalidN.value || invalidK.value || invalidPerms.value)
const invalidMsg = computed(() => [
  invalidN.value ? t('web.analysis.consensus_bad_n') : '',
  invalidK.value ? t('web.analysis.consensus_bad_k') : '',
  invalidPerms.value ? t('web.analysis.consensus_bad_perms') : ''
].filter(Boolean).join('；'))

const consJob = useJob((res) => {
  const r = res as unknown as ModulesConsensusResult
  result.value = r
  if (r.ok) emit('saved', r.files_written ?? [])
})

function run(): void {
  if (consJob.running.value) return
  if (anyInvalid.value) return  // 拒绝提交:参数行下方的红字实时提示就是反馈
  // 空框回写默认值(可见化,不静默)
  nRef.value = effN.value
  kRef.value = effK.value
  permsRef.value = effPerms.value
  result.value = null
  void consJob.run(
    api.modulesConsensus({
      directory_rel: props.ppi.directory_rel,
      prefix: props.ppi.prefix,
      mode: props.mode,
      drug: props.drug,
      disease: props.disease,
      genes_text: props.genesText,
      n: effN.value,
      k: effK.value,
      n_perms: effPerms.value
    }),
    { doneText: t('web.analysis.consensus_done') }
  )
}

/* 模板窄化别名(vue-tsc 对嵌套 v-if 的窄化不可靠,同 AnalysisView reachCmp 模式) */
const cons = computed(() => {
  const r = result.value
  return r && r.ok ? r : null
})
const failReason = computed(() => {
  const r = result.value
  return r && !r.ok ? r.reason ?? '' : ''
})
const reasonText = computed(() => {
  if (!failReason.value) return ''
  const ik = `web.analysis.consensus_reason_${failReason.value}`
  const s = t(ik)
  return s === ik ? failReason.value : s
})

const seedLine = computed(() =>
  props.mode === 'names'
    ? t('web.analysis.consensus_seed_disease', { disease: props.disease })
    : t('web.analysis.consensus_seed_genes')
)

const shown = computed(() => (cons.value?.rows ?? []).slice(0, 100))
const fillGenes = computed(() => cons.value?.consistent_list ?? [])

/* ---- viridis 顺序色系(锚点色与 output/modules_panel_demo 静态样例同源) ---- */
const VIRIDIS = ['#440154', '#46327e', '#365c8d', '#277f8e', '#1fa187', '#4ac16d', '#a0da39', '#fde725']
function hexRgb(h: string): [number, number, number] {
  return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)]
}
function virRgb(tt: number): [number, number, number] {
  const x = Math.min(Math.max(tt, 0), 1) * (VIRIDIS.length - 1)
  const i = Math.min(Math.floor(x), VIRIDIS.length - 2)
  const f = x - i
  const a = hexRgb(VIRIDIS[i])
  const b = hexRgb(VIRIDIS[i + 1])
  return [Math.round(a[0] + (b[0] - a[0]) * f), Math.round(a[1] + (b[1] - a[1]) * f), Math.round(a[2] + (b[2] - a[2]) * f)]
}
function virStyle(rank: number | null, nTotal: number): Record<string, string> {
  if (rank === null || rank < 1) return {}
  const rgb = virRgb((rank - 1) / Math.max(nTotal - 1, 1))
  const lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
  return {
    background: `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`,
    color: lum < 140 ? 'var(--c-on-accent)' : 'var(--c-ink)'
  }
}
const legendSwatches = computed(() => [0, 1 / 3, 2 / 3, 1].map((tt) => {
  const rgb = virRgb(tt)
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`
}))
function rankTitle(rank: number | null, nTotal: number): string {
  return rank === null
    ? t('web.analysis.consensus_not_member')
    : t('web.analysis.consensus_rank_title', { r: rank, n: nTotal })
}
</script>

<template>
  <section class="sec cons-panel">
    <h3 class="sec-title">{{ t('web.analysis.consensus_title') }}</h3>
    <p class="hint-line">{{ t('web.analysis.consensus_hint') }}</p>
    <p class="hint-line">{{ seedLine }}</p>

    <div class="cons-params">
      <label class="cons-field">N
        <input v-model.number="nRef" type="number" min="10" max="300" :class="{ bad: invalidN }" />
      </label>
      <label class="cons-field">k
        <input v-model.number="kRef" type="number" min="1" max="3" :class="{ bad: invalidK }" />
      </label>
      <label class="cons-field">n_perms
        <input v-model.number="permsRef" type="number" min="0" max="200" :class="{ bad: invalidPerms }" />
      </label>
      <button class="cons-run" :disabled="consJob.running.value" @click="run">
        {{ consJob.running.value ? t('web.analysis.rerunning') : t('web.analysis.consensus_run') }}
      </button>
    </div>
    <p v-if="anyInvalid" class="notice err">{{ invalidMsg }}</p>
    <!-- 运行进度(治"死态"):细线进度条与主分析同款全局样式;后端步进上报
         (DIAMOnD 每入模一步 / 显著性每完成一次置换)驱动连续爬行 -->
    <div v-if="consJob.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((consJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
    </div>
    <p v-if="consJob.running.value" class="hint-line">
      {{ consJob.text.value }}<template v-if="consJob.fraction.value !== null"> · {{ Math.round(consJob.fraction.value * 100) }}%</template>
    </p>
    <p v-else-if="consJob.error.value" class="notice err">{{ consJob.error.value }}</p>

    <p v-if="failReason" class="notice err">{{ reasonText }}</p>
    <template v-if="cons">
      <div class="cons-chips">
        <span class="chip">{{ t('web.analysis.consensus_chip_seeds', { i: cons.coverage.seeds_in_lcc, total: cons.coverage.seeds_total }) }}</span>
        <span class="chip">{{ t('web.analysis.consensus_chip_sizes', { s: `${cons.module_sizes.rwr}/${cons.module_sizes.diamond}/${cons.module_sizes.neighbor}` }) }}</span>
        <span class="chip">{{ t('web.analysis.consensus_chip_consistent', { c: cons.consistent_genes, total: cons.rows_total }) }}</span>
        <span class="chip">{{ t('web.analysis.consensus_chip_k', { k: cons.params.k, n: cons.params.n }) }}</span>
        <span v-if="!cons.significance" class="chip dim">{{ t('web.analysis.consensus_chip_nosig') }}</span>
      </div>

      <template v-if="cons.significance">
        <p class="hint-line">{{ t('web.analysis.consensus_sig_hint') }}</p>
        <div class="tbl-wrap">
          <table class="tbl">
            <thead><tr>
              <th>{{ t('web.analysis.consensus_col_method') }}</th>
              <th>p</th>
              <th>n_perms</th>
              <th>{{ t('web.analysis.consensus_sig_obs') }}</th>
              <th>{{ t('web.analysis.consensus_sig_null') }}</th>
            </tr></thead>
            <tbody>
              <tr v-for="(s, m) in cons.significance" :key="m">
                <td>{{ m }}</td>
                <td>{{ fmt(s.p, 4) }}</td>
                <td>{{ s.n_perms }}</td>
                <td>{{ s.observed_cohesion }}</td>
                <td>{{ fmt(s.null_mean, 1) }} ± {{ fmt(s.null_sd, 1) }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </template>

      <div class="tbl-wrap">
        <table class="tbl cons-tbl">
          <thead><tr>
            <th>#</th>
            <th>{{ t('web.analysis.consensus_col_gene') }}</th>
            <th class="r">RWR</th>
            <th class="r">DIAMOnD</th>
            <th class="r">{{ t('web.analysis.consensus_col_neighbor') }}</th>
            <th class="r">{{ t('web.analysis.consensus_col_score') }}</th>
            <th class="r">{{ t('web.analysis.consensus_col_methods') }}</th>
            <th>{{ t('web.analysis.consensus_col_badge') }}</th>
          </tr></thead>
          <tbody>
            <tr v-for="(row, i) in shown" :key="row.gene">
              <td class="ix">{{ i + 1 }}</td>
              <td class="gn">{{ row.gene }}</td>
              <td class="rk" :class="{ none: row.rwr_rank === null }" :style="virStyle(row.rwr_rank, cons.params.n)" :title="rankTitle(row.rwr_rank, cons.params.n)">{{ row.rwr_rank ?? '–' }}</td>
              <td class="rk" :class="{ none: row.diamond_rank === null }" :style="virStyle(row.diamond_rank, cons.params.n)" :title="rankTitle(row.diamond_rank, cons.params.n)">{{ row.diamond_rank ?? '–' }}</td>
              <td class="rk" :class="{ none: row.neighbor_rank === null }" :style="virStyle(row.neighbor_rank, cons.params.n)" :title="rankTitle(row.neighbor_rank, cons.params.n)">{{ row.neighbor_rank ?? '–' }}</td>
              <td class="sc">{{ fmt(row.consensus_score, 1) }}</td>
              <td class="nm">{{ row.n_methods }}</td>
              <td>
                <span v-if="row.is_seed" class="bd seed">{{ t('web.analysis.consensus_badge_seed') }}</span>
                <span v-if="row.consistent" class="bd cons">{{ t('web.analysis.consensus_badge_consistent') }}</span>
                <span v-if="!row.is_seed && !row.consistent" class="bd solo">{{ t('web.analysis.consensus_badge_single') }}</span>
              </td>
            </tr>
          </tbody>
        </table>
      </div>

      <p class="hint-line lg">
        {{ t('web.analysis.consensus_legend') }}
        <span v-for="(sw, i) in legendSwatches" :key="i" class="sw" :style="{ background: sw }" />
        <span class="none-sw">–</span>
      </p>
      <div class="cons-foot">
        <span class="hint-line">{{ t('web.analysis.consensus_rows_hint', { total: cons.rows_total, shown: shown.length }) }}</span>
        <!-- 退化态(种子 ≥ N,三方法未扩张):一致子集=种子本身,填入是无意义的原样
             回倒——禁用按钮并明示原因(§1.12 陷阱②) -->
        <button v-if="!cons.degenerate" class="cons-fill" @click="emit('fill', fillGenes)">
          {{ t('web.analysis.consensus_fill_btn', { count: fillGenes.length }) }}
        </button>
        <span v-else class="hint-line blocked">{{ t('web.analysis.consensus_fill_blocked', { s: cons.coverage.seeds_in_lcc, n: cons.params.n }) }}</span>
      </div>
    </template>
  </section>
</template>

<style scoped>
/* chrome 色值全走 var(--c-*)(深浅主题自适应);viridis 色格是数据编码,
 * 深浅主题同一套(见文件头注释) */
.cons-panel {
  margin-top: 34px;
}
.cons-params {
  display: flex;
  align-items: center;
  gap: 14px;
  margin: 10px 0 4px;
  flex-wrap: wrap;
}
.cons-field {
  display: inline-flex;
  align-items: baseline;
  gap: 6px;
  font-size: 12.5px;
  color: var(--c-ink2);
}
.cons-field input {
  width: 64px;
  font: inherit;
  font-size: 13px;
  color: var(--c-ink);
  background: var(--c-panel);
  border: none;
  border-bottom: 1px solid var(--c-line);
  padding: 2px 4px;
}
.cons-field input:focus {
  outline: none;
  border-bottom-color: var(--c-accent);
}
.cons-field input.bad {
  border-bottom-color: var(--c-danger);
}
.cons-run,
.cons-fill {
  background: var(--c-accent);
  color: var(--c-on-accent);
  border: none;
  font: inherit;
  font-size: 13px;
  padding: 7px 16px;
  border-radius: 3px;
  cursor: pointer;
}
.blocked {
  max-width: 60%;
  text-align: right;
}
.cons-run:disabled {
  opacity: 0.5;
  cursor: default;
}
.cons-chips {
  display: flex;
  gap: 8px;
  margin: 12px 0 4px;
  flex-wrap: wrap;
}
.chip {
  font-size: 12px;
  color: var(--c-ink);
  background: var(--c-rowline);
  border: 1px solid var(--c-line);
  padding: 3px 10px;
  border-radius: 999px;
}
.chip.dim {
  color: var(--c-mid);
}
.cons-tbl th.r,
.cons-tbl td.rk,
.cons-tbl td.sc,
.cons-tbl td.nm {
  text-align: center;
}
.cons-tbl .ix {
  color: var(--c-ink3);
  width: 30px;
}
.cons-tbl .gn {
  font-weight: 600;
  color: var(--c-ink);
  width: 110px;
}
.cons-tbl .rk {
  width: 84px;
  font-variant-numeric: tabular-nums;
}
.cons-tbl .rk.none {
  color: var(--c-ink4);
  background: var(--c-hover);
}
.cons-tbl .sc {
  color: var(--c-mid2);
  font-variant-numeric: tabular-nums;
  width: 70px;
}
.cons-tbl .nm {
  color: var(--c-body);
  width: 60px;
}
.bd {
  display: inline-block;
  font-size: 11px;
  padding: 1px 8px;
  border-radius: 999px;
}
.bd.cons {
  background: var(--c-accent);
  color: var(--c-on-accent);
}
.bd.seed {
  border: 1px solid var(--c-accent);
  color: var(--c-ink);
  margin-right: 4px;
}
.bd.solo {
  border: 1px solid var(--c-flag-off);
  color: var(--c-ink2);
}
.lg .sw {
  display: inline-block;
  width: 26px;
  height: 12px;
  margin-left: 6px;
}
.lg .none-sw {
  color: var(--c-ink4);
  margin-left: 6px;
}
.cons-foot {
  display: flex;
  justify-content: space-between;
  align-items: center;
  border-top: 1px solid var(--c-line);
  margin-top: 14px;
  padding-top: 12px;
  gap: 12px;
  flex-wrap: wrap;
}
</style>
