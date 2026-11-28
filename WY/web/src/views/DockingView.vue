<script setup lang="ts">
/** 分子对接(阶段 3b,Molecular_Docking 的 TS 迁移):六步向导。
 * 会话状态存服务端(core.docking_session),前端持 sid;网络重步骤走 job 轮询;
 * 3D 预览为服务端 py3Dmol 渲染的 HTML,iframe 内嵌。
 * Phase 4 拆分后壳层只管:会话与 sid、Step 2 配体、Step 5 运行、跨步共享状态
 * (s1/receptor/box/poseSel/dockState)、两个 reset 编排器与 busyLock 聚合;
 * Step 1/3/4/6 各为 DockingSourcePanel/ReceptorPanel/BoxPanel/ResultPanel,
 * 经 v-model(receptor/box/poseSel)+ 事件 + expose(clear/busy/gene/桥接位)接线。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type DockingSourceSummary, type DockingStructCheck, type Row } from '../api'
import DockingResultPanel from '../components/DockingResultPanel.vue'
import DockingBoxPanel from '../components/DockingBoxPanel.vue'
import DockingReceptorPanel from '../components/DockingReceptorPanel.vue'
import DockingSourcePanel from '../components/DockingSourcePanel.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'

const { t } = useI18n()

/* ---------------- 会话 ---------------- */

const sid = ref('')
/* Step5 参数行的实时值;拉取失败时用内置默认兜底 */
const dockParams = ref<{ e: number; m: number; c: number; timeout: number }>(
  { e: 8, m: 9, c: 4, timeout: 600 })

async function ensureSession(): Promise<string> {
  if (sid.value) return sid.value
  sid.value = (await api.dockingSession()).sid
  return sid.value
}

/* ---------------- Step 1 数据源与候选基因(组件 DockingSourcePanel) ----------------
 * s1 与结构筛选结果留壳层:genes 计算、Step 2 药名语境与 Step 6 归档前后缀都要
 * 读;三来源切换/结构筛选随组件。组件 emit picked/cleared/checked 回传,壳层
 * 负责 s1 写入、下游清场与配体自动获取(Step 2 状态在壳层)。 */

const s1 = ref<DockingSourceSummary | null>(null)
const structResult = ref<DockingStructCheck | null>(null)
const srcPanel = ref<InstanceType<typeof DockingSourcePanel> | null>(null)
/* 受体候选过滤开关:勾选框在来源面板,请求在受体面板,壳层桥接 */
const allowAfBridge = computed(() => srcPanel.value?.allowAf ?? false)
// 手动会话:Step 2 隐藏"按药名查 PubChem"(无药名语境;CID 直查与 SMILES 保留)
const isManualSession = computed(() => !!srcPanel.value?.manualSession)

const genes = computed<string[]>(() => structResult.value?.genes ?? s1.value?.candidates ?? [])

/** 来源就绪(标准来源选择 autoLigand=true / 手动基因提交 false):写 s1、清下游;
 * 标准路径再自动获取配体(镜像 Streamlit:识别到药物名后自动 PubChem 一次,
 * 失败可手动改后重试;手动路径不动已填的查询串)。 */
function onSourcePicked(p: DockingSourceSummary & { autoLigand: boolean }): void {
  s1.value = { drug: p.drug, disease: p.disease, hub: p.hub, enr: p.enr, candidates: p.candidates }
  structResult.value = null
  resetDownstream()
  if (p.autoLigand) {
    ligQuery.value = p.drug || ''
    if (p.drug) void fetchLigand()
  }
}

/** 来源三态切换:清当前来源与下游,避免两个流程共用半截会话 */
function onSourceCleared(): void {
  s1.value = null
  structResult.value = null
  resetDownstream()
}

/** 结构筛选完成(job onDone 或"不做筛选"的直接路径) */
function onStructChecked(r: DockingStructCheck): void {
  structResult.value = r
}

/** 数据源变更后清空前端本地的下游视图状态(服务端会话已由端点同步失效)。
 * 组件化后各步组件 expose clear(),壳层只编排顺序、不再直接摸子组件状态。 */
function resetDownstream(): void {
  ligStale.value = false
  lig.value = null
  receptor.value = null
  dockState.value = null
  poseSel.value = 1
  recPanel.value?.clear()
  boxPanel.value?.clear()
  resultPanel.value?.clear()
}

/* ---------------- Step 2 配体 ---------------- */

const ligSource = ref<'pubchem' | 'smiles'>('pubchem')
const ligQuery = ref('')
const ligSmiles = ref('')
const lig = ref<{ ok: boolean; error?: string; cid?: number | null; lig_pdbqt?: string } | null>(null)
const ligBusy = ref(false)
const ligStale = ref(false)

async function fetchLigand(): Promise<void> {
  sid.value = await ensureSession()
  // 代际守卫(P0-B):发起时捕获来源上下文快照,期间切换则丢弃,防旧配体写回
  // 新上下文(来源状态在 DockingSourcePanel 内,经 expose 的 ligCtx 取签名)
  const ctx0 = srcPanel.value?.ligCtx ?? ''
  ligStale.value = false
  ligBusy.value = true
  try {
    const r = ligSource.value === 'pubchem'
      ? await api.dockingLigand(sid.value, { source: 'pubchem', query: ligQuery.value })
      : await api.dockingLigand(sid.value, { source: 'smiles', smiles: ligSmiles.value })
    if ((srcPanel.value?.ligCtx ?? '') !== ctx0) {
      ligStale.value = true
      return
    }
    lig.value = r
  } catch (e) {
    lig.value = { ok: false, error: e instanceof Error ? e.message : String(e) }
  } finally {
    ligBusy.value = false
  }
}

/* ---------------- Step 3 受体(组件 DockingReceptorPanel) ----------------
 * receptor 壳层持有:canRun、Step 3 完成徽标与对接盒面板都要读(v-model 双向);
 * 检索/上传两条入路、候选表与制备 job 随组件。组件 emit reset-downstream 时
 * 由下方 resetSteps4to6 编排清场。 */

const receptor = ref<{ pdb_id: string; receptor_source: string } | null>(null)
const recPanel = ref<InstanceType<typeof DockingReceptorPanel> | null>(null)
/* 归档面板的结论占位与 3D 图注需要当前基因(状态在受体面板内,经 expose 桥接) */
const recGene = computed(() => recPanel.value?.gene ?? '')

/** 换基因 / 换受体后清空 Step 4-6 的旧结果:盒子、对接得分、归档、3D 视图
 * 都绑定在具体受体上,留着上一受体的缓存展示会误导。 */
function resetSteps4to6(): void {
  dockState.value = null
  poseSel.value = 1
  boxPanel.value?.clear()
  resultPanel.value?.clear()
}

/* ---------------- Step 4 对接盒(组件 DockingBoxPanel) ----------------
 * box 壳层持有:canRun 与 Step 4 完成徽标要读;推断/微调/预览随组件走(v-model)。 */

const box = ref<{ center: number[]; size: number[] } | null>(null)
const boxPanel = ref<InstanceType<typeof DockingBoxPanel> | null>(null)

/* ---------------- Step 5 运行对接 ---------------- */

const runSeq = ref(0)  // 每完成一轮对接 +1,归档签名用它识别"新一轮结果未归档"
const runJob = useJob((result) => {
  const r = result as { log_tail: string; scores: Row[] }
  dockState.value = { result: { log_tail: r.log_tail, scores: r.scores } }
  runSeq.value++
})
const dockState = ref<{ result: { log_tail: string; scores: Row[] } } | null>(null)

const canRun = computed(() => !!(
  lig.value?.ok && receptor.value && box.value))

// 任一步骤 job 或配体获取进行中,锁住来源/基因切换入口(P0-B):
// 中途换来源/基因会让在途结果落到新上下文,拼出「新来源 + 旧配体/旧受体」组合
// (结构筛选/受体制备的 busy 分属来源/受体面板,经 expose 聚合)
const busyLock = computed(() =>
  ligBusy.value || runJob.running.value ||
  !!srcPanel.value?.busy || !!recPanel.value?.busy)

async function runDock(): Promise<void> {
  sid.value = await ensureSession()
  dockState.value = null
  await runJob.run(api.dockingRun(sid.value), { doneText: t('web.docking.run_done_text') })
}

/* ---------------- Step 6 结果与归档(组件 DockingResultPanel) ----------------
 * poseSel 壳层持有:Step 5 得分表行点击与归档面板构象下拉共写一个值(v-model)。
 * 归档签名/3D/互作全部随组件;受体 pdb_id 与 runSeq 经 props 喂给它。 */

const poseSel = ref(1)
const resultPanel = ref<InstanceType<typeof DockingResultPanel> | null>(null)

const scores = computed<Row[]>(() => dockState.value?.result?.scores ?? [])
const logTail = computed(() => dockState.value?.result?.log_tail ?? '')

/* ---------------- 通用 ---------------- */

const { notice, flash } = useNotice(3600)

/* 完成徽标只余壳层直管的 Step 2/5;Step 1/3/4/6 随各自组件的自身状态计算 */
const stepDone = computed(() => ({
  s2: !!lig.value?.ok,
  s5: scores.value.length > 0
}))

onMounted(async () => {
  // Step5 参数行取实时值(exhaustiveness/cpu 面板可调,超时随深度放宽);失败退回
  // 内置默认(来源清单拉取随 DockingSourcePanel 自己的 onMounted)
  try {
    const s = await api.settings()
    dockParams.value = { e: s.exhaustiveness, m: s.num_modes, c: s.docking_cpu,
                         timeout: s.docking_timeout }
  } catch {
    /* 保持内置默认显示 */
  }
})
</script>

<template>
  <p class="kicker s1">{{ t('web.docking.kicker') }}</p>
  <h1 class="title s2">{{ t('web.docking.title') }}</h1>
  <p class="sub s3">{{ t('web.docking.sub') }}</p>

  <div class="panel wide s4">
    <div class="xlinks">
      <a class="bar-btn" href="#/stars">{{ t('web.nav.stars_link') }}</a>
    </div>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <!-- ============ Step 1(组件 DockingSourcePanel:s1/结构结果留壳层,picked/cleared/checked 回传) ============ -->
    <DockingSourcePanel ref="srcPanel" :busy-lock="busyLock" :s1="s1" :struct-result="structResult"
                        :ensure-session="ensureSession" :flash="flash"
                        @picked="onSourcePicked" @cleared="onSourceCleared" @checked="onStructChecked" />

    <!-- ============ Step 2 ============ -->
    <section class="step" :class="{ done: stepDone.s2 }">
      <h3 class="step-title">{{ t('web.docking.s2_title') }} <span v-if="stepDone.s2" class="ok-mark">{{ t('web.docking.s2_ready') }}{{ lig?.cid ? ` · CID ${lig.cid}` : ' · SMILES' }}</span></h3>
      <div class="lv-group">
        <button class="lv" :class="{ on: ligSource === 'pubchem' }" @click="ligSource = 'pubchem'">PubChem</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: ligSource === 'smiles' }" @click="ligSource = 'smiles'">SMILES</button>
      </div>
      <template v-if="ligSource === 'pubchem'">
        <div class="set-row">
          <label class="set-label">{{ isManualSession ? t('web.docking.lig_cid_manual') : t('web.docking.lig_label') }}</label>
          <div class="field"><input v-model="ligQuery" :placeholder="isManualSession ? t('web.docking.cid_ph') : (s1?.drug || t('web.docking.drug_name_ph'))" spellcheck="false" /><i class="line"></i></div>
          <button class="bar-btn" :disabled="ligBusy" @click="fetchLigand">
            {{ ligBusy ? t('web.docking.fetching') : t('web.docking.fetch_ligand') }}
          </button>
        </div>
        <p v-if="!isManualSession && !s1?.drug && !lig?.ok" class="hint-line">{{ t('web.docking.no_drug_hint') }}</p>
      </template>
      <template v-else>
        <textarea v-model="ligSmiles" rows="2" class="ta" :placeholder="t('web.docking.smiles_ph')"></textarea>
        <button class="bar-btn" :disabled="ligBusy || !ligSmiles.trim()" @click="fetchLigand">
          {{ ligBusy ? t('web.docking.generating') : t('web.docking.gen_conformer') }}
        </button>
      </template>
      <p v-if="ligStale" class="stale">{{ t('web.job.stale') }}</p>
      <p v-if="lig?.ok" class="hint-line good">{{ t('web.docking.lig_ok') }}</p>
      <p v-if="lig && !lig.ok" class="stale">{{ lig.error }}</p>
    </section>

    <!-- ============ Step 3(组件 DockingReceptorPanel:receptor 双向,下游失效走 emit) ============ -->
    <DockingReceptorPanel ref="recPanel" v-model:receptor="receptor" :genes="genes"
                          :allow-af="allowAfBridge" :ensure-session="ensureSession" :flash="flash"
                          @reset-downstream="resetSteps4to6" />

    <!-- ============ Step 4(组件 DockingBoxPanel:box 双向,canRun 需要) ============ -->
    <DockingBoxPanel ref="boxPanel" v-model:box="box" :sid="sid" :receptor="receptor"
                     :ensure-session="ensureSession" :flash="flash" />

    <!-- ============ Step 5 ============ -->
    <section class="step" :class="{ done: stepDone.s5 }">
      <h3 class="step-title">{{ t('web.docking.s5_title') }} <span v-if="stepDone.s5" class="ok-mark">✓</span></h3>
      <p v-if="canRun" class="hint-line">{{ t('web.docking.run_params', dockParams) }}</p>
      <button v-if="canRun" class="bar-btn strong" :disabled="runJob.running.value" @click="runDock">
        {{ runJob.running.value ? t('web.docking.running_vina') : t('web.docking.run_btn') }}
      </button>
      <p v-else class="hint-line">{{ t('web.docking.need_s234') }}</p>
      <div v-if="runJob.running.value" class="progress">
        <div class="progress-track">
          <div class="progress-fill indeterminate"></div>
        </div>
        <p class="progress-text">{{ runJob.text.value }}</p>
      </div>
      <p v-if="runJob.error.value" class="stale">{{ t('web.docking.run_fail', { msg: runJob.error.value }) }}</p>
      <template v-if="scores.length">
        <p class="hint-line">{{ t('web.docking.log_tail_hint') }}</p>
        <pre class="log">{{ logTail }}</pre>
        <div class="tbl-wrap">
          <table class="tbl">
            <thead><tr><th v-for="k in Object.keys(scores[0])" :key="k">{{ k }}</th></tr></thead>
            <tbody>
              <tr v-for="(r, i) in scores" :key="i" :class="{ selrow: poseSel === Number(r.mode) }"
                  @click="poseSel = Number(r.mode)">
                <td v-for="(v, k) in r" :key="k">{{ v }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </template>
    </section>

    <!-- ============ Step 6(组件 DockingResultPanel:poseSel 双向) ============ -->
    <DockingResultPanel ref="resultPanel" v-model:poseSel="poseSel" :sid="sid" :scores="scores"
                        :pdb-id="receptor?.pdb_id ?? ''" :gene="recGene" :drug="s1?.drug ?? ''"
                        :disease="s1?.disease ?? ''" :run-seq="runSeq"
                        :ensure-session="ensureSession" :flash="flash" />
  </div>
</template>

