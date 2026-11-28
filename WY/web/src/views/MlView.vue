<script setup lang="ts">
/** GEO 表达谱 + ML 共识筛选(阶段 3c,Expression_ML_Screening 的 TS 迁移)。
 * 忠实迁移原页面四区布局:① 输入(GEO 检索/直填/本地矩阵) ② 配置(分组正则/排序表/参数)
 * ③ 运行(七阶段 job) ④ 结果(指标/ROC/共识/漏斗/shortlist + 近 3 次摘要)。
 * 末段并入 ML_Results_Information(结果管理):run 清单/搜索/删除/ZIP/keep-N 清理。
 * 全链编排唯一定义在 core.ml_pipeline.run_screen_chain(CLI/页面/API 共用)。 */
import { computed, onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type MlHistoryRun, type MlOverviewResult,
         type MlProbeResult, type MlRunResult, type Row } from '../api'
import MlCountsPicker from '../components/MlCountsPicker.vue'
import MlGroupPicker from '../components/MlGroupPicker.vue'
import MlResultPanel from '../components/MlResultPanel.vue'
import Stepper from '../components/Stepper.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import { dropBadge, pruneBadges } from '../utils/probeBadgeCache'
import { useMlGseProbe } from '../composables/useMlGseProbe'
import { useMlProbeBadges } from '../composables/useMlProbeBadges'
import { useMlRun } from '../composables/useMlRun'
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
      resetProbe()
      // 徽章缓存不清:判定不可变,换检索词同 GSE 即时复用;只剪过期条目
      pruneBadges()
      selectedRows.value = new Set()
  } catch (e) {
    flash(t('web.ml.search_fail', { msg: e instanceof Error ? e.message : String(e) }))
  } finally {
    searchBusy.value = false
  }
}

/** 点击检索表行 = 选中该 GSE(与对接页受体选择同一交互):驱动预检面板与运行配置,
 * 行内 GSE 编号即确认信息;已选行再点保持选中(无反选)。 */
function pickRow(g: string): void {
  if (!g) return
  if (gse.value.trim().toUpperCase() === g) return
  gse.value = g
  onGseChange()
}

const searchStale = computed(() =>
  searchResults.value.length > 0 &&
  disease.value.trim().toLowerCase() !== searchDisease.value.toLowerCase())

/* ---------------- ② 配置 ---------------- */

const casePat = ref('')
const ctrlPat = ref('')
const overviewJob = useJob((result) => {
  // 徽章不作废(2026-11-28 定稿):总览只读数据集,归档推断在 TTL 内仍可信,
  // 跑一次总览就作废过于激进;作废时机移至 runJob onDone(真正落盘数据集)
  overview.value = result as unknown as MlOverviewResult
  flash(t('web.ml.overview_done_flash'))
})
const overview = ref<MlOverviewResult | null>(null)

/* counts 列名全集:供 MlCountsPicker 算 spill(建议前缀波及组外列数) */
const countsSampleCols = computed(() =>
  overview.value?.counts_rows?.map((r) => String(r.sample ?? '')) ?? [])

/* 分组选择器填入的公共落点(转义后的正则由组件 emit,壳层只写两个输入框) */
function applyCaseCtrl(caseRe: string, ctrlRe: string): void {
  casePat.value = caseRe
  ctrlPat.value = ctrlRe
}

async function loadOverview(): Promise<void> {
  const g = gse.value.trim().toUpperCase()
  if (!g) { flash(t('web.ml.need_gse')); return }
  // 代际守卫(P0-F):发起时捕获 GSE 快照,期间切换则总览作废(否则旧数据集的
  // group_picker 回填 ② 区,gpApply/ctsApply 会把旧分组正则填进新 GSE 的 run)
  await overviewJob.run(api.mlOverview(g), {
    doneText: t('web.ml.overview_done_text'),
    stale: () => gse.value.trim().toUpperCase() !== g
  })
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

const forceJob = useJob(() => {
  // 失效钩子:强制重下已落盘,旧徽章立即作废(gse 取值须在清 fallback 之前)
  dropBadge(fallback.value?.gse ?? gse.value.trim().toUpperCase())
  fallback.value = null
  resetProbe()
  flash(t('web.ml.force_done_flash'))
})

async function forceRedownload(): Promise<void> {
  const g = fallback.value?.gse ?? gse.value.trim().toUpperCase()
  if (!g) return
  await forceJob.run(api.mlForceRefresh(g), { doneText: t('web.ml.force_done_text') })
}

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

/* ① 区逻辑组合:徽章批检(勾选/圆点/批检)+ GSE 单检与可用性锁;
 * searchPlatformFor 由批检组产出、单检消费。onGseChange/pickRow/search/
 * forceRedownload 经 resetProbe 重置单检态。 */
const { selectedRows, probeBatchBusy, selectedCount, allSelected, rowGse, toggleRow,
        toggleAll, symbolVerdict, dotClass, searchPlatformFor, probeBatch } =
  useMlProbeBadges(searchResults, flash)
const { exprPath, labelsPath, localCase, probe, probeBusy, badgeApplied, probeGse,
        onGseChange, resetProbe, useLocal, localPartial, hasGse, sourceConflict,
        probeR, probeErr, probeAssess, unusable, fmtEta } =
  useMlGseProbe(gse, searchPlatformFor)

/* ③ 区运行编排(run 提交/七阶段显示/在途恢复):参数包与守卫经 deps 注入;
 * 强制重下(forceJob)跨①区失效钩子(resetProbe),留壳层。 */
const { resultData, fallback, runJob, canRun, run, resumeRunJob,
        STAGE_ORDER, STATE_ICON, stageLabels, stageStates } = useMlRun({
  gse, disease, tab, casePat, ctrlPat, overview,
  exprPath, labelsPath, localCase, rankedSel, topK, prefilterN, bootstrapN,
  hasGse, useLocal, localPartial, sourceConflict, unusable,
  flash, loadHistory
})

onMounted(() => {
  void loadDiseases()
  void loadRanked()
  void loadHistory()
  void resumeRunJob()
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
          <div v-if="searchResults.length" class="probe-row">
            <span class="hint-line">{{ t('web.ml.probe_sel_hint', { n: selectedCount }) }}</span>
            <button class="bar-btn" :disabled="probeBatchBusy || !selectedCount"
                    @click="probeBatch($event.shiftKey)">
              {{ probeBatchBusy ? t('web.ml.probe_sel_running') : t('web.ml.probe_sel_btn') }}
            </button>
            <span class="hint-line">{{ t('web.ml.probe_force_hint') }}</span>
          </div>
          <div v-if="searchResults.length" class="tbl-wrap">
            <p v-if="searchStale" class="stale">{{ t('web.ml.search_stale', { d: searchDisease }) }}</p>
            <table class="tbl">
              <thead><tr>
                <th class="col-chk"><input type="checkbox" :checked="allSelected" @change="toggleAll" /></th>
                <th v-for="k in Object.keys(searchResults[0])" :key="k">{{ k }}</th>
              </tr></thead>
              <tbody>
                <tr v-for="(r, i) in searchResults" :key="i"
                    :class="{ selrow: rowGse(r) && rowGse(r) === gse.trim().toUpperCase() }"
                    style="cursor: pointer" @click="pickRow(rowGse(r))">
                  <td class="col-chk" @click.stop><input type="checkbox" :checked="selectedRows.has(rowGse(r))"
                                                        @change="toggleRow(rowGse(r))" /></td>
                  <td v-for="(v, k) in r" :key="k">
                    <template v-if="k === 'gse'">
                      <span class="dot" :class="dotClass(rowGse(r))"
                            :title="t('web.ml.probe_dot_' + dotClass(rowGse(r)))"></span>{{ v }}
                    </template>
                    <template v-else-if="k === 'symbol_level'">
                      <span v-if="symbolVerdict(rowGse(r))"
                            :class="symbolVerdict(rowGse(r))!.tone"
                            :title="symbolVerdict(rowGse(r))!.tone === 'unknown' ? t('web.ml.sym_unknown_title') : ''">
                        {{ t(symbolVerdict(rowGse(r))!.key) }}
                      </span>
                      <span v-else class="sym-hint">{{ v }}</span>
                    </template>
                    <template v-else>{{ v }}</template>
                  </td>
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
        <div class="field">
          <input v-model="gse" placeholder="GSE63067" spellcheck="false"
                 @change="onGseChange" /><i class="line"></i>
        </div>
        <span v-if="gse" class="state on">{{ t('web.ml.confirmed', { gse }) }}</span>
        <span class="hint-line">{{ t('web.ml.pick_row_hint') }}</span>
      </div>

      <!-- 预检展示(徽章豁免时显示占位说明,2026-11-10 QA) -->
      <template v-if="tab !== 'local' && gse && (probeBusy || probe || badgeApplied)">
        <p v-if="probeBusy" class="hint-line">{{ t('web.ml.probing') }}</p>
        <p v-else-if="badgeApplied" :class="badgeApplied.status === 'unusable' ? 'stale' : 'hint-line good'">
          {{ badgeApplied.status === 'unusable'
              ? t('web.ml.badge_unusable', { n: badgeApplied.ageMin })
              : t('web.ml.badge_skip_probe', { n: badgeApplied.ageMin }) }}
        </p>
        <template v-else>
          <p v-if="probeR" class="hint-line">
            {{ t('web.ml.probe_meta', { size: fmtSize(probeR.size_bytes),
                                        speed: (probeR.speed_bps / 1024).toFixed(0),
                                        eta: fmtEta(Math.min(probeR.size_bytes / Math.max(probeR.speed_bps, 1), 5400)) }) }}
          </p>
          <p v-if="probeR && ((probeR.size_bytes >= 1048576 && !probeR.reliable) || probeR.size_bytes / Math.max(probeR.speed_bps, 1) > 300)"
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
        <!-- probeBusy 置灰:可用性判定出炉前不启动分钟级下载,免得对不可用数据集白跑(QA 2026-11-09) -->
        <button class="bar-btn" :disabled="!gse || overviewJob.running.value || unusable || probeBusy"
                @click="loadOverview">
          {{ overviewJob.running.value ? t('web.ml.overview_running') : t('web.ml.overview_btn') }}
        </button>
        <span v-if="overviewJob.stale.value" class="state">{{ t('web.job.stale') }}</span>
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
        <!-- counts 解析失败独立通道:后端正常完成的 job 也可能带 counts_error
             (counts 文件损坏/格式异常),不渲染的话用户会把失败误判成"没有 counts" -->
        <p v-if="overview.counts_error" class="stale">
          {{ t('web.ml.counts_error', { msg: overview.counts_error }) }}
        </p>
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
          <!-- counts 模式分组选择器:同外壳不同发动机(见 core.expression.counts_picker_fields) -->
          <p v-if="overview.counts_picker && !overview.counts_picker.fields?.length" class="stale">
            {{ t('web.ml.gp_counts_nomatch',
                 { m: overview.counts_match?.matched ?? 0, n: overview.counts_match?.total ?? 0 }) }}
          </p>
          <template v-if="overview.counts_picker?.fields?.length">
            <MlCountsPicker :fields="overview.counts_picker.fields" :all-cols="countsSampleCols"
                            :gse="gse" :unusable="unusable" :flash="flash" @apply="applyCaseCtrl" />
          </template>
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
          <template v-if="overview.group_picker?.fields?.length">
            <MlGroupPicker :fields="overview.group_picker.fields"
                           :gse="gse" :unusable="unusable" :flash="flash" @apply="applyCaseCtrl" />
          </template>
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
      <!-- probeBusy 置灰是防御性双保险:换 GSE 会清空分组正则,canRun 理论上已为假;
           显式挡住防止未来动线变化后在可用性判定出炉前启动分钟级 run(QA 2026-11-09) -->
      <button class="build" :disabled="!canRun || runJob.running.value || overviewJob.running.value || probeBusy"
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
      <!-- 过期提示必须在 running 块外:useJob 的过期文案写在 text 里,随进度块一起消失 -->
      <p v-if="runJob.stale.value" class="stale">{{ t('web.job.stale') }}</p>

      <div v-if="fallback" class="set-row" style="margin-top:14px">
        <span class="stale">{{ t('web.ml.fallback_warn', { date: fallback.date }) }}</span>
        <button class="bar-btn" :disabled="forceJob.running.value" @click="forceRedownload">
          {{ forceJob.running.value ? t('web.ml.forcing') : t('web.ml.force_btn') }}
        </button>
      </div>
      <p v-if="forceJob.error.value" class="stale">{{ t('web.ml.force_fail', { msg: forceJob.error.value }) }}</p>
    </section>

    <!-- ============ ④ 结果区 ============ -->
    <MlResultPanel :result="resultData" :flash="flash" />

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
  color: var(--c-ink);
  font-weight: 500;
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
.build:disabled {
  opacity: 0.4;
  cursor: default;
}
.build-hint {
  font-size: 11px;
  letter-spacing: 0.3em;
  color: var(--c-ink3);
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
  color: var(--c-ink4);
}
.stage.running {
  color: var(--c-ink);
}
.stage.done {
  color: var(--c-ink2);
}
.stage-icon {
  margin-right: 6px;
  font-style: normal;
}
.stale-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.2em;
  color: var(--c-danger);
  border: 1px solid var(--c-danger);
  padding: 0 5px;
  margin-left: 10px;
  vertical-align: 1.5px;
}
/* 平台排序表:定高滚动列表(发丝线边框,一行一项,候选多时滚动选) */
.ranked-list {
  margin-top: 10px;
  max-height: 176px;
  overflow-y: auto;
  border: 1px solid var(--c-line-soft);
}
.ranked-item {
  display: flex;
  align-items: baseline;
  gap: 10px;
  padding: 7px 12px;
  border-bottom: 1px solid var(--c-rowline);
  font-size: 12px;
  letter-spacing: 0.04em;
  color: var(--c-body);
  cursor: pointer;
  user-select: none;
}
.ranked-item:last-child {
  border-bottom: none;
}
.ranked-item:hover {
  background: var(--c-hover);
}
.ranked-item input {
  accent-color: var(--c-ink);
  flex: none;
}
.ranked-name {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-variant-numeric: tabular-nums;
}
/* 排序表列表与下方参数调节器之间留出呼吸间距;
   三个参数组用 flex 拉开互相间隔(列间距 64px),窄屏自动换行(行距 20px) */
.opt-grid {
  margin-top: 30px;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 20px 64px;
}
/* 每组内"标签文字 + 步进器"绑成一个整体,文字与 +/- 按钮留足间距 */
.opt {
  display: inline-flex;
  align-items: center;
  gap: 16px;
}
.opt-label {
  white-space: nowrap;
}
</style>
