<script setup lang="ts">
/** ML ④ 结果区:指标卡 / ROC / 共识·漏斗·入围三表 + 下载区(整包 zip +
 * 产物归类列表)。纯展示组件,全部数据从 result(MlRunResult)投影派生;
 * 整包下载自持,失败提示经 props.flash 走页顶 notice(父级 useNotice 实例)。
 * 产物归类:页面展示对应的四项为主成果;DEG 全表/共识曲线/元数据等为补充
 * 产物(默认折叠)。未识别的新产物文件名即说明,进补充区——后端新增产物
 * 不需要改前端就能出现。 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl, type MlRunResult, type Row } from '../api'

const props = defineProps<{
  result: MlRunResult | null
  flash: (msg: string) => void
}>()
const { t } = useI18n()

const meta = computed(() => props.result?.meta ?? null)
const consensusRows = computed<Row[]>(() => props.result?.consensus ?? [])
const consensusTotal = computed(() => props.result?.consensus_total ?? 0)
const funnelRows = computed<Row[]>(() => props.result?.funnel ?? [])
const shortlistRows = computed<Row[]>(() => props.result?.shortlist ?? [])
const shortlistCount = computed(() => props.result?.shortlist_count ?? '—')
const artifactsMap = computed<Record<string, string>>(() => props.result?.artifacts ?? {})

/* 整包下载:当前 run 目录打成一个 zip(单 run 走 make_dir_zip) */
const zipBusy = ref(false)
const runDirName = computed(() => (props.result?.out_dir ?? '').split(/[\\/]/).pop() ?? '')
async function zipRun(): Promise<void> {
  if (!props.result?.out_dir) return
  zipBusy.value = true
  try {
    const blob = await api.mlRunsZip([props.result.out_dir])
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${runDirName.value}.zip`
    a.click()
    URL.revokeObjectURL(url)
  } catch {
    props.flash(t('web.common.zip_fail'))
  } finally {
    zipBusy.value = false
  }
}

interface DlItem { file: string; rel: string; desc: string; main: boolean }

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
</script>

<template>
  <section class="sec" v-if="result">
    <h3 class="sec-title">{{ t('web.ml.sec_result', { gse: result.gse, s: result.elapsed }) }}</h3>
    <div class="metrics">
      <div class="metric"><b>{{ meta?.n_samples }}</b><span>{{ t('web.ml.m_samples') }}</span></div>
      <div class="metric"><b>{{ meta?.n_genes }}</b><span>{{ t('web.ml.m_genes') }}</span></div>
      <div class="metric">
        <b>{{ meta?.cv_mean ?? meta?.cv_auc ?? 'N/A' }}</b>
        <span>CV-AUC<template v-if="meta?.cv_sd != null"> ± {{ meta.cv_sd }}</template></span>
      </div>
      <div class="metric"><b>{{ meta?.algo_spearman_mean }}</b><span>{{ t('web.ml.m_spearman') }}</span></div>
    </div>

    <template v-if="result.roc_png_rel">
      <p class="hint-line" style="margin-top:16px">{{ t('web.ml.roc_hint') }}</p>
      <p v-if="(meta?.n_samples ?? 0) < 10" class="hint-line">{{ t('web.ml.roc_small_n') }}</p>
      <img class="fig-img" :src="artifactsUrl(result.roc_png_rel)" alt="ROC" />
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
</template>

<style scoped>
.fig-img {
  max-width: 420px;
  margin-top: 10px;
  border: 1px solid var(--c-line-soft);
  background: var(--c-panel);
}
.sub-sec {
  margin: 22px 0 8px;
  font-size: 13px;
  font-weight: 500;
  letter-spacing: 0.12em;
  color: var(--c-ink);
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
  color: var(--c-quiet);
}
.fold[open] .fold2::before {
  content: '⌄ ';
}
</style>
