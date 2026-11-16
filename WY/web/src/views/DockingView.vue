<script setup lang="ts">
/** 分子对接(阶段 3b,Molecular_Docking 的 TS 迁移):六步向导。
 * 会话状态存服务端(core.docking_session),前端持 sid;网络重步骤走 job 轮询;
 * 3D 预览为服务端 py3Dmol 渲染的 HTML,iframe 内嵌。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl, type Row } from '../api'
import Stepper from '../components/Stepper.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'

const { t } = useI18n()

/* ---------------- 会话 ---------------- */

const sid = ref('')
const dockError = ref('')
/* Step5 参数行的实时值;拉取失败时用内置默认兜底 */
const dockParams = ref<{ e: number; m: number; c: number; timeout: number }>(
  { e: 8, m: 9, c: 4, timeout: 600 })

async function ensureSession(): Promise<string> {
  if (sid.value) return sid.value
  sid.value = (await api.dockingSession()).sid
  return sid.value
}

/* ---------------- Step 1 数据源与候选基因 ---------------- */

const srcMode = ref<'ppi' | 'enr'>('ppi')
const ppiDirs = ref<{ name: string; display: string }[]>([])
const enrFiles = ref<{ name: string; display: string }[]>([])
const selSource = ref('')
const manualGenes = ref('')
// 手动基因对接专属入口:manual_* 桶不再混入 PPI/富集下拉,由此处接管
const manualMode = ref(false)
const manualDirs = ref<{ name: string; n_genes: number; mtime: string; genes: string[] }[]>([])
const manualPick = ref('')
const s1 = ref<{ drug: string; disease: string; hub: number; enr: number; candidates: string[] } | null>(null)
// 手动会话:Step 2 隐藏"按药名查 PubChem"(无药名语境;CID 直查与 SMILES 保留)
const isManualSession = computed(() => manualMode.value && !s1.value?.drug)
const useStruct = ref(true)
const allowAf = ref(false)
const structJob = useJob((result) => {
  structResult.value = result as { genes: string[]; no_struct: string[]; hub_no_struct: string[]; structures: Record<string, number> }
})
const structResult = ref<{ genes: string[]; no_struct: string[]; hub_no_struct: string[]; structures: Record<string, number> } | null>(null)

const genes = computed<string[]>(() => structResult.value?.genes ?? s1.value?.candidates ?? [])

async function loadSources(): Promise<void> {
  const s = await api.dockingSources()
  ppiDirs.value = s.ppi_dirs
  enrFiles.value = s.enr_files
  manualDirs.value = s.manual_dirs ?? []
}
async function selectSource(name: string): Promise<void> {
  sid.value = await ensureSession()
  const r = await api.dockingSource(sid.value, { mode: srcMode.value, name })
  s1.value = { drug: r.drug, disease: r.disease, hub: r.hub_count, enr: r.enr_count, candidates: r.candidates }
  structResult.value = null
  resetDownstream()
  // 自动获取配体(镜像 Streamlit:识别到药物名后自动 PubChem 一次,失败可手动改后重试)
  ligQuery.value = r.drug || ''
  if (r.drug) void fetchLigand()
}

function pickManualFile(name: string): void {
  const d = manualDirs.value.find((x) => x.name === name)
  if (d?.genes?.length) {
    manualGenes.value = d.genes.join(', ')
    flash(t('web.docking.manual_loaded', { n: d.genes.length }))
  }
}

/** 手动面板输入框 → 基因列表(与后端同款分隔符切分,实时计数) */
const parsedManualGenes = computed<string[]>(() =>
  manualGenes.value.split(/[,;\s]+/).map((g) => g.trim().toUpperCase()).filter(Boolean))

/** 来源三态切换(ppi/enr/manual):入口级操作,清掉当前来源与下游,
 * 避免两个流程共用半截会话;首次切到手动时补拉清单。 */
function setSrcMode(mode: 'ppi' | 'enr' | 'manual'): void {
  if (mode === 'manual') {
    manualMode.value = true
    s1.value = null
    structResult.value = null
    resetDownstream()
    if (!manualDirs.value.length) void loadSources()
    return
  }
  manualMode.value = false
  if (srcMode.value === mode) return
  srcMode.value = mode
  s1.value = null
  structResult.value = null
  resetDownstream()
}

async function useManual(): Promise<void> {
  if (!manualGenes.value.trim()) return
  sid.value = await ensureSession()
  const r = await api.dockingSource(sid.value, { mode: 'manual', manual_genes: manualGenes.value })
  s1.value = { drug: r.drug, disease: r.disease, hub: r.hub_count, enr: r.enr_count, candidates: r.candidates }
  structResult.value = null
  resetDownstream()
}

/** 数据源变更后清空前端本地的下游视图状态(服务端会话已由端点同步失效)。 */
function resetDownstream(): void {
  ligStale.value = false
  lig.value = null
  receptor.value = null
  candRows.value = []
  candError.value = ''
  pdbSel.value = ''
  gene.value = ''
  pickedFile.value = null
  uploadReport.value = null
  box.value = null
  boxMeta.value = null
  dockState.value = null
  saved.value = null
  starred.value = false
  poseSel.value = 1
  conclusion.value = ''
  view3dSrc.value = ''
}

async function checkStructures(): Promise<void> {
  sid.value = await ensureSession()
  const body = { use_struct: useStruct.value, allow_af: allowAf.value }
  if (!useStruct.value) {
    // 不做结构筛选:候选直接作为最终靶点,无 job
    const r = await api.dockingCheckStructures(sid.value, body)
    structResult.value = (r.direct ?? {}) as typeof structResult.value
    return
  }
  await structJob.run(api.dockingCheckStructures(sid.value, body) as Promise<{ job_id: string }>,
    { doneText: t('web.docking.check_done_text') })
}

const candidatePreview = computed(() => {
  const list = s1.value?.candidates ?? []
  return list.slice(0, 60).join(', ') + (list.length > 60 ? ' …' : '')
})

/* ---------------- Step 2 配体 ---------------- */

const ligSource = ref<'pubchem' | 'smiles'>('pubchem')
const ligQuery = ref('')
const ligSmiles = ref('')
const lig = ref<{ ok: boolean; error?: string; cid?: number | null; lig_pdbqt?: string } | null>(null)
const ligBusy = ref(false)
const ligStale = ref(false)

async function fetchLigand(): Promise<void> {
  sid.value = await ensureSession()
  // 代际守卫(P0-B):发起时捕获来源上下文快照,期间切换则丢弃,防旧配体写回新上下文
  const ctx0 = `${srcMode.value}|${selSource.value}|${manualMode.value ? 'manual' : 'auto'}`
  ligStale.value = false
  ligBusy.value = true
  try {
    const r = ligSource.value === 'pubchem'
      ? await api.dockingLigand(sid.value, { source: 'pubchem', query: ligQuery.value })
      : await api.dockingLigand(sid.value, { source: 'smiles', smiles: ligSmiles.value })
    if (`${srcMode.value}|${selSource.value}|${manualMode.value ? 'manual' : 'auto'}` !== ctx0) {
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

/* ---------------- Step 3 受体 ---------------- */

const gene = ref('')
const candRows = ref<{ pdb_id: string; resolution: string; method: string; title: string; source: string; flag: string; note: string }[]>([])
const candError = ref('')
const pdbSel = ref('')
const receptor = ref<{ pdb_id: string; receptor_source: string } | null>(null)
/** 受体来源:RCSB 检索 / 自带结构上传(两条入路,下游盒子-Vina-归档共用) */
const recSource = ref<'rcsb' | 'upload'>('rcsb')
const pickedFile = ref<File | null>(null)
const uploadReport = ref<{ n_atom: number; n_protein: number; n_nucleic: number;
                           n_other_atom: number; n_lig_het: number; n_water: number;
                           chains: string[]; sha256: string; size_bytes: number;
                           warnings: string[] } | null>(null)

const recJob = useJob((result) => {
  const r = result as unknown as { pdb_id: string; receptor_source: string;
                                   report?: typeof uploadReport.value }
  receptor.value = r
  uploadReport.value = r.report ?? null
  resetSteps4to6()
  flash(t(r.report ? 'web.docking.upload_ready' : 'web.docking.receptor_ready'))
})

const recFailLabel = computed(() =>
  recSource.value === 'upload' ? t('web.docking.upload_fail') : t('web.docking.prep_fail'))

/** 切换受体来源:清 receptor 与下游(两条入路的结果互斥),候选表保留便于切回 */
function setRecSource(mode: 'rcsb' | 'upload'): void {
  if (recSource.value === mode || recJob.running.value) return
  recSource.value = mode
  receptor.value = null
  pdbSel.value = ''
  pickedFile.value = null
  uploadReport.value = null
  resetSteps4to6()
}

function onFilePicked(e: Event): void {
  const input = e.target as HTMLInputElement
  const f = input.files?.[0] ?? null
  input.value = '' // 允许重复选择同一个文件重新上传
  if (!f) return
  pickedFile.value = f
  uploadReport.value = null
  if (receptor.value) {  // 换文件即换结构:旧受体作废
    receptor.value = null
    resetSteps4to6()
  }
}

async function uploadReceptor(): Promise<void> {
  if (!gene.value || !pickedFile.value) return
  sid.value = await ensureSession()
  // 代际守卫(P0-B):基因切换或切回 RCSB 分支都视为上下文失效
  const g0 = gene.value
  await recJob.run(api.dockingReceptorUpload(sid.value, gene.value, pickedFile.value),
    { doneText: t('web.docking.upload_ready'),
      stale: () => gene.value !== g0 || recSource.value !== 'upload' })
}

/** 换基因 / 换受体后清空 Step 4-6 的旧结果:盒子、对接得分、归档、3D 视图
 * 都绑定在具体受体上,留着上一受体的缓存展示会误导。 */
function resetSteps4to6(): void {
  box.value = null
  boxMeta.value = null
  dockState.value = null
  saved.value = null
  starred.value = false
  poseSel.value = 1
  conclusion.value = ''
  view3dSrc.value = ''
  view3dKey.value++
  interRows.value = []
  pocketRows.value = []
  interOpen.value = false
}

async function loadCandidates(): Promise<void> {
  if (!gene.value) return
  sid.value = await ensureSession()
  pdbSel.value = ''
  receptor.value = null
  resetSteps4to6()
  try {
    const r = await api.dockingReceptorCandidates(sid.value, { gene: gene.value, allow_af: allowAf.value })
    candError.value = r.error ?? ''
    candRows.value = r.rows
  } catch (e) {
    candError.value = e instanceof Error ? e.message : String(e)
    candRows.value = []
  }
}

function flagBlocked(flag: string): boolean {
  return flag === 'dna_only'
}

async function prepareReceptor(): Promise<void> {
  if (!pdbSel.value || flagBlocked(selFlag.value)) return
  sid.value = await ensureSession()
  // 代际守卫(P0-B):发起时捕获基因快照,期间切换则受体结果作废,
  // 不落 Step 3 展示(服务端另有写前核对兜底)
  const g0 = gene.value
  await recJob.run(api.dockingReceptorPrepare(sid.value, { gene: gene.value, pdb_id: pdbSel.value }),
    { doneText: t('web.docking.receptor_ready'),
      stale: () => gene.value !== g0 })
}

const selFlag = computed(() => candRows.value.find((r) => r.pdb_id === pdbSel.value)?.flag ?? '')

/* ---------------- Step 4 对接盒 ---------------- */

const box = ref<{ center: number[]; size: number[] } | null>(null)
const boxMeta = ref<{ clamped: boolean; inferred_ligands: string[] | null; residues: string[] } | null>(null)
const boxBusy = ref(false)
// 盒子预览重建计数:iframe 的 src 固定指向 box-3d,微调/重推不重建就永远不刷新
// (QA 2026-11-09 反馈)。递增同时改 src(?v=)与 :key,防缓存且强制重挂载。
const boxKey = ref(0)

async function inferBox(): Promise<void> {
  sid.value = await ensureSession()
  boxBusy.value = true
  try {
    const r = await api.dockingBox(sid.value, { infer: true })
    box.value = r.box
    boxMeta.value = { clamped: r.clamped, inferred_ligands: r.inferred_ligands, residues: r.residues }
    boxKey.value++
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  } finally {
    boxBusy.value = false
  }
}

async function applyManualBox(): Promise<void> {
  sid.value = await ensureSession()
  boxBusy.value = true
  try {
    const r = await api.dockingBox(sid.value, { infer: false, center: box.value!.center, size: box.value!.size })
    boxMeta.value = { clamped: false, inferred_ligands: null, residues: r.residues }
    boxKey.value++
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  } finally {
    boxBusy.value = false
  }
}

/** 盒坐标输入解析(2026-11 审计 M2):空串 Number('')===0 会把坐标静默清零,
 * '-'/'.' 等逐字符中间态是 NaN 会被后端 422 拒掉——两类都直接放弃本次编辑,
 * 保留旧值,等用户输完整数字再落。 */
function parseBoxInput(e: Event): number | null {
  const raw = (e.target as HTMLInputElement).value
  if (raw.trim() === '') return null
  const v = Number(raw)
  return Number.isFinite(v) ? v : null
}

function editCenter(i: number, e: Event): void {
  const v = parseBoxInput(e)
  if (v === null || !box.value) return
  box.value = { center: box.value.center.map((c, j) => (j === i ? v : c)),
                size: box.value.size }
}
function editSize(i: number, e: Event): void {
  const v = parseBoxInput(e)
  if (v === null || !box.value) return
  box.value = { center: box.value.center,
                size: box.value.size.map((s, j) => (j === i ? v : s)) }
}

const fmt = (v: unknown, digits = 3): string =>
  typeof v === 'number' ? v.toFixed(digits) : (v == null ? '—' : String(v))

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
const busyLock = computed(() =>
  ligBusy.value || structJob.running.value || recJob.running.value || runJob.running.value)

async function runDock(): Promise<void> {
  sid.value = await ensureSession()
  dockState.value = null
  await runJob.run(api.dockingRun(sid.value), { doneText: t('web.docking.run_done_text') })
}

/* ---------------- Step 6 结果与归档 ---------------- */

const poseSel = ref(1)
const conclusion = ref('')
const drugName = ref('')
const diseaseName = ref('')
const saved = ref<{ result_dir: string } | null>(null)
const starred = ref(false)
const saveBusy = ref(false)

const scores = computed<Row[]>(() => dockState.value?.result?.scores ?? [])
const logTail = computed(() => dockState.value?.result?.log_tail ?? '')

/* 归档防重复:签名 = 对接轮次 + 受体 + 构象 + 结论 + 前后缀;归档后按钮置灰,
 * 上述任一内容变化才重新允许归档(想留多版本就改一下结论再归)。 */
const lastSavedSig = ref('')
const archiveSig = computed(() => JSON.stringify(
  [runSeq.value, receptor.value?.pdb_id, poseSel.value,
   conclusion.value, drugName.value, diseaseName.value]))
const archiveDirty = computed(() => archiveSig.value !== lastSavedSig.value)

const cutoff = ref(4.0)
const focus = ref(false)
const view3dKey = ref(0)
const view3dSrc = ref('')
const view3dBusy = ref(false)
const pocketNames = computed(() => pocketRows.value.map((r) => `${r.resn}${r.resi}`))
const pocketRows = ref<{ resn: string; resi: string; chain: string; dist: number }[]>([])
const interRows = ref<Row[]>([])
const interOpen = ref(false)

async function load3d(): Promise<void> {
  if (!scores.value.length || !sid.value) return
  view3dBusy.value = true
  try {
    view3dKey.value++  // 重建 iframe:同参数再点"刷新视图"也能真正重载,而不是静默无反应
    view3dSrc.value = api.dockingView3dUrl(sid.value, poseSel.value, cutoff.value, focus.value)
    const inter = await api.dockingInteractions(sid.value, poseSel.value)
    interRows.value = inter.rows
    pocketRows.value = inter.rows.map((r) => ({
      resn: String(r.resn), resi: String(r.resi), chain: String(r.chain ?? ''),
      dist: Number(r.dist ?? 0)
    }))
  } catch {
    view3dSrc.value = ''
  } finally {
    view3dBusy.value = false
  }
}

let reloadTimer = 0
function scheduleReload(): void {
  window.clearTimeout(reloadTimer)
  reloadTimer = window.setTimeout(() => void load3d(), 300)
}

async function saveResult(): Promise<void> {
  sid.value = await ensureSession()
  saveBusy.value = true
  try {
    const r = await api.dockingSave(sid.value, {
      pose_mode: poseSel.value, conclusion: conclusion.value,
      drug: drugName.value, disease: diseaseName.value
    })
    saved.value = { result_dir: r.result_dir }
    lastSavedSig.value = archiveSig.value
    starred.value = false
    flash(t('web.docking.archived'))
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  } finally {
    saveBusy.value = false
  }
}

async function toggleStar(): Promise<void> {
  if (!saved.value) return
  try {
    starred.value = (await api.dockingStar(sid.value, !starred.value)).starred
  } catch {
    flash(t('web.docking.star_fail'))
  }
}

const save3dBusy = ref(false)
async function save3d(): Promise<void> {
  if (!saved.value) return
  save3dBusy.value = true
  try {
    const caption = `NPAP · ${drugName.value || '-'} × ${diseaseName.value || '-'} · ${gene.value} · pose ${poseSel.value}`
    const r = await api.dockingSave3d(sid.value,
      { pose: poseSel.value, caption, cutoff: cutoff.value })
    flash(t('web.docking.saved_3d', { rel: r.artifact_rel,
                                       offline: r.offline_ok ? t('web.docking.offline_ok') : t('web.docking.offline_no') }))
    // 保存进结果目录的同时把 HTML 直接下载到本地(对齐旧版"下载交互式 3D 视图")
    const a = document.createElement('a')
    a.href = artifactsUrl(r.artifact_rel, true)
    a.click()
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  } finally {
    save3dBusy.value = false
  }
}

/* ---------------- 通用 ---------------- */

const { notice, flash } = useNotice(3600)

const stepDone = computed(() => ({
  s1: genes.value.length > 0,
  s2: !!lig.value?.ok,
  s3: !!receptor.value,
  s4: !!box.value,
  s5: scores.value.length > 0,
  s6: !!saved.value
}))

onMounted(async () => {
  try {
    await loadSources()
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  }
  // Step5 参数行取实时值(exhaustiveness/cpu 面板可调,超时随深度放宽);失败退回内置默认
  try {
    const s = await api.settings()
    dockParams.value = { e: s.exhaustiveness, m: s.num_modes, c: s.docking_cpu,
                         timeout: s.docking_timeout }
  } catch {
    /* 保持内置默认显示 */
  }
})
onBeforeUnmount(() => {
  window.clearTimeout(reloadTimer)
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

    <!-- ============ Step 1 ============ -->
    <section class="step" :class="{ done: stepDone.s1 }">
      <h3 class="step-title">{{ t('web.docking.s1_title') }} <span v-if="stepDone.s1" class="ok-mark">{{ t('web.docking.s1_genes_ok', { n: genes.length }) }}</span></h3>

      <!-- 三来源切换:PPI / 富集缓存 / 手动基因(专属面板,下游 Step 2-6 共用) -->
      <div class="lv-group">
        <button class="lv" :class="{ on: srcMode === 'ppi' }" @click="setSrcMode('ppi')">{{ t('web.docking.src_ppi') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: srcMode === 'enr' }" @click="setSrcMode('enr')">{{ t('web.docking.src_enr') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: manualMode }" @click="setSrcMode('manual')">{{ t('web.docking.manual_entry_btn') }}</button>
      </div>

      <!-- 手动基因面板 -->
      <template v-if="manualMode">
        <p class="hint-line">{{ t('web.docking.manual_panel_hint') }}</p>
        <div class="set-row">
          <label class="set-label">{{ t('web.docking.manual_pick_label') }}</label>
          <select class="sel" v-model="manualPick" @change="pickManualFile(manualPick)"
                  :disabled="!manualDirs.length">
            <option value="" disabled>{{ t('web.docking.manual_pick_ph') }}</option>
            <option v-for="d in manualDirs" :key="d.name" :value="d.name">
              {{ d.n_genes }} {{ t('web.docking.manual_gene_unit') }} · {{ d.mtime }} · {{ d.name }}
            </option>
          </select>
        </div>
        <textarea v-model="manualGenes" rows="3" class="ta"
                  :placeholder="t('web.docking.manual_genes_ph')" spellcheck="false"></textarea>
        <p class="hint-line">{{ t('web.docking.manual_count', { n: parsedManualGenes.length }) }}</p>
        <button class="bar-btn" :disabled="!parsedManualGenes.length || structJob.running.value"
                @click="useManual">{{ t('web.docking.manual_btn') }}</button>
      </template>

      <!-- 标准流程来源选择 -->
      <template v-else>
      <p v-if="srcMode === 'ppi' && !ppiDirs.length" class="hint-line">{{ t('web.docking.no_ppi_dirs') }}</p>
      <select v-if="srcMode === 'ppi' && ppiDirs.length" class="sel" v-model="selSource"
              :disabled="busyLock" @change="selectSource(selSource)">
        <option value="" disabled>{{ t('web.docking.sel_dir_ph') }}</option>
        <option v-for="d in ppiDirs" :key="d.name" :value="d.name">{{ d.display }}</option>
      </select>
      <select v-if="srcMode === 'enr' && enrFiles.length" class="sel" v-model="selSource"
              :disabled="busyLock" @change="selectSource(selSource)">
        <option value="" disabled>{{ t('web.docking.sel_enr_ph') }}</option>
        <option v-for="f in enrFiles" :key="f.name" :value="f.name">{{ f.display }}</option>
      </select>
      </template>

      <template v-if="s1">
        <div class="metrics three">
          <div v-if="!isManualSession" class="metric"><b>{{ s1.hub }}</b><span>{{ t('web.docking.m_hub') }}</span></div>
          <div v-if="!isManualSession" class="metric"><b>{{ s1.enr }}</b><span>{{ t('web.docking.m_enr') }}</span></div>
          <div class="metric"><b>{{ s1.candidates.length }}</b><span>{{ t('web.docking.m_cand') }}</span></div>
        </div>
        <p class="hint-line">{{ candidatePreview }}</p>

        <div class="filters">
          <label class="chk"><input v-model="useStruct" type="checkbox" :disabled="!s1.candidates.length" /> {{ t('web.docking.chk_struct') }}</label>
          <label class="chk"><input v-model="allowAf" type="checkbox" /> {{ t('web.docking.chk_af') }}</label>
          <button v-if="useStruct && s1.candidates.length" class="bar-btn"
                  :disabled="structJob.running.value" @click="checkStructures">
            {{ structJob.running.value ? t('web.docking.checking') : t('web.docking.check_struct_btn') }}
          </button>
        </div>
        <div v-if="structJob.running.value" class="progress">
          <div class="progress-track">
            <div class="progress-fill" :style="{ width: `${((structJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
          </div>
          <p class="progress-text">{{ structJob.text.value }}</p>
        </div>
        <p v-if="structJob.error.value" class="stale">{{ t('web.docking.check_fail', { msg: structJob.error.value }) }}</p>
        <template v-if="structResult">
          <p class="hint-line good">
            {{ t('web.docking.has_struct', { n: structResult.genes.length,
                                             list: structResult.genes.slice(0, 40).join(', ') }) }}
          </p>
          <p v-if="structResult.no_struct.length" class="hint-line">
            {{ t('web.docking.no_struct_removed', { n: structResult.no_struct.length,
                                                    list: structResult.no_struct.slice(0, 40).join(', ') }) }}
          </p>
          <p v-if="structResult.hub_no_struct.length" class="stale">
            {{ t('web.docking.hub_removed_warn', { n: structResult.hub_no_struct.length,
                                                   names: structResult.hub_no_struct.join(', ') }) }}
          </p>
        </template>
      </template>
      <p v-if="!manualMode && !s1 && !ppiDirs.length && !enrFiles.length" class="hint-line">
        {{ t('web.docking.no_sources_hint') }}
      </p>
    </section>

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

    <!-- ============ Step 3 ============ -->
    <section class="step" :class="{ done: stepDone.s3 }">
      <h3 class="step-title">{{ t('web.docking.s3_title') }} <span v-if="stepDone.s3" class="ok-mark">✓ {{ receptor?.pdb_id }}</span></h3>
      <template v-if="genes.length">
        <select class="sel" v-model="gene" :disabled="recJob.running.value" @change="loadCandidates">
          <option value="" disabled>{{ t('web.docking.sel_gene_ph') }}</option>
          <option v-for="g in genes" :key="g" :value="g">{{ g }}</option>
        </select>
        <!-- 受体来源:RCSB 检索 / 自带结构上传(同源模建、MD 快照、突变体等) -->
        <div class="lv-group">
          <button class="lv" :class="{ on: recSource === 'rcsb' }" :disabled="recJob.running.value"
                  @click="setRecSource('rcsb')">{{ t('web.docking.rec_src_rcsb') }}</button>
          <i class="lv-sep">/</i>
          <button class="lv" :class="{ on: recSource === 'upload' }" :disabled="recJob.running.value"
                  @click="setRecSource('upload')">{{ t('web.docking.rec_src_upload') }}</button>
        </div>
        <template v-if="recSource === 'rcsb'">
          <p v-if="candError" class="stale">{{ candError }}</p>
          <template v-if="candRows.length">
            <div class="tbl-wrap">
              <table class="tbl">
                <thead><tr><th>PDB</th><th>{{ t('web.docking.th_res') }}</th><th>{{ t('web.docking.th_method') }}</th><th>{{ t('web.docking.th_title') }}</th><th>{{ t('web.docking.th_note') }}</th></tr></thead>
                <tbody>
                  <tr v-for="r in candRows" :key="r.pdb_id" :class="{ selrow: pdbSel === r.pdb_id }"
                      @click="!flagBlocked(r.flag) && (pdbSel = r.pdb_id)">
                    <td>{{ r.pdb_id }}</td><td>{{ r.resolution }}</td><td>{{ r.method }}</td>
                    <td class="ttl">{{ r.title }}</td>
                    <td :class="{ warn: !!r.note }">{{ r.note }}</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p v-if="selFlag === 'dna_only'" class="stale">{{ t('web.docking.dna_only') }}</p>
            <button class="bar-btn" :disabled="!pdbSel || flagBlocked(selFlag) || recJob.running.value"
                    @click="prepareReceptor">
              {{ recJob.running.value ? t('web.docking.prepping') : t('web.docking.prep_btn') }}
            </button>
          </template>
          <p v-else-if="gene && !candError && !recJob.running.value" class="hint-line">
            {{ t('web.docking.rcsb_empty') }}
          </p>
        </template>
        <template v-else>
          <p class="hint-line">{{ t('web.docking.upload_hint') }}</p>
          <div class="upload-row">
            <label class="bar-btn upload-btn">
              <input type="file" accept=".pdb,.ent,.cif,.mmcif" @change="onFilePicked" />
              {{ t('web.docking.upload_pick') }}
            </label>
            <span v-if="pickedFile" class="hint-line">{{ pickedFile.name }}</span>
          </div>
          <button class="bar-btn" :disabled="!gene || !pickedFile || recJob.running.value"
                  @click="uploadReceptor">
            {{ recJob.running.value ? t('web.docking.uploading') : t('web.docking.upload_btn') }}
          </button>
        </template>
        <!-- 进度/错误/stale 为两条入路共用 -->
        <div v-if="recJob.running.value" class="progress">
          <div class="progress-track">
            <div class="progress-fill" :style="{ width: `${((recJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
          </div>
          <p class="progress-text">{{ recJob.text.value }}</p>
        </div>
        <p v-if="recJob.error.value" class="notice err">{{ recFailLabel }}:{{ recJob.error.value }}</p>
        <p v-if="recJob.stale.value" class="stale">{{ t('web.job.stale') }}</p>
        <!-- 上传体检报告:自动清洗说明 + 逐条隐患提示(后端中文文案原样展示) -->
        <template v-if="uploadReport">
          <p class="hint-line">
            {{ t('web.docking.upload_report', { n: uploadReport.n_protein + uploadReport.n_nucleic + uploadReport.n_other_atom,
                                                lig: uploadReport.n_lig_het, water: uploadReport.n_water }) }}
          </p>
          <p v-for="(w, i) in uploadReport.warnings" :key="i" class="stale">{{ w }}</p>
        </template>
        <p v-if="receptor?.receptor_source === 'alphafold'" class="stale">
          {{ t('web.docking.af_warn', { pdb: receptor.pdb_id.split(':')[1] }) }}
        </p>
      </template>
      <p v-else class="hint-line">{{ t('web.docking.need_s1') }}</p>
    </section>

    <!-- ============ Step 4 ============ -->
    <section class="step" :class="{ done: stepDone.s4 }">
      <h3 class="step-title">{{ t('web.docking.s4_title') }} <span v-if="stepDone.s4" class="ok-mark">✓</span></h3>
      <button class="bar-btn" :disabled="!receptor || boxBusy" @click="inferBox">
        {{ boxBusy ? t('web.docking.inferring') : t('web.docking.infer_box') }}
      </button>
      <p v-if="boxMeta?.clamped && receptor?.receptor_source === 'upload'" class="stale">
        {{ t('web.docking.upload_no_lig_note') }}
      </p>
      <p v-else-if="boxMeta?.clamped" class="stale">{{ t('web.docking.clamped_note') }}</p>
      <p v-else-if="boxMeta?.inferred_ligands?.length" class="hint-line">
        {{ t('web.docking.ref_ligands', { names: boxMeta.inferred_ligands.join(', ') }) }}
      </p>
      <template v-if="box">
        <div class="box-grid">
          <div class="box-nums">
            <p class="ctl-label">{{ t('web.docking.box_label') }}</p>
            <div class="num-grid">
              <label v-for="(c, i) in box.center" :key="'c' + i" class="num">
                <span>c{{ 'xyz'[i] }}</span>
                <input type="number" step="0.5" :value="c" @change="editCenter(i, $event)" />
              </label>
              <label v-for="(s, i) in box.size" :key="'s' + i" class="num">
                <span>s{{ 'xyz'[i] }}</span>
                <input type="number" step="1" :value="s" @change="editSize(i, $event)" />
              </label>
            </div>
            <button class="bar-btn" :disabled="boxBusy" @click="applyManualBox">{{ t('web.docking.apply_tune') }}</button>
            <p v-if="boxMeta?.residues?.length" class="hint-line">
              {{ t('web.docking.box_residues', { n: boxMeta.residues.length,
                                                  list: boxMeta.residues.slice(0, 30).join(', ') }) }}
            </p>
            <p v-else class="hint-line">{{ t('web.docking.box_empty') }}</p>
          </div>
          <div class="box-view">
            <iframe :key="boxKey" class="fig-frame" :src="`/api/docking/${sid}/box-3d?v=${boxKey}`" title="box preview" />
          </div>
        </div>
      </template>
    </section>

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

    <!-- ============ Step 6 ============ -->
    <section class="step" :class="{ done: stepDone.s6 }">
      <h3 class="step-title">{{ t('web.docking.s6_title') }}
        <span v-if="stepDone.s6" class="ok-mark">✓ {{ saved?.result_dir.split('/').pop() || saved?.result_dir }}</span>
      </h3>
      <template v-if="scores.length">
        <div class="set-row">
          <label class="set-label">{{ t('web.docking.pose_label') }}</label>
          <select class="sel inline" v-model.number="poseSel" @change="scheduleReload">
            <option v-for="r in scores" :key="String(r.mode)" :value="Number(r.mode)">mode {{ r.mode }}</option>
          </select>
          <span class="ctl-meta">{{ t('web.docking.affinity', { v: scores.find((r) => Number(r.mode) === poseSel)?.affinity }) }}</span>
        </div>
        <div class="set-row">
          <label class="set-label">{{ t('web.docking.conclusion_label') }}</label>
          <textarea v-model="conclusion" rows="2" class="ta wide"
                    :placeholder="t('web.docking.conclusion_ph', { gene, aff: scores[0]?.affinity })"></textarea>
        </div>
        <div class="set-row">
          <label class="set-label">{{ t('web.docking.prefix_label') }}</label>
          <div class="field"><input v-model="drugName" :placeholder="s1?.drug || t('web.docking.drug_name_ph')" spellcheck="false" /><i class="line"></i></div>
          <div class="field"><input v-model="diseaseName" :placeholder="s1?.disease || t('web.docking.disease_name_ph')" spellcheck="false" /><i class="line"></i></div>
          <button class="bar-btn strong" :disabled="saveBusy || (!!saved && !archiveDirty)" @click="saveResult">
            {{ saveBusy ? t('web.docking.archiving')
              : (!!saved && !archiveDirty) ? t('web.docking.archived_btn') : t('web.docking.archive_btn') }}
          </button>
          <button v-if="saved" class="bar-btn" @click="toggleStar">
            {{ starred ? t('web.docking.starred') : t('web.docking.star') }}
          </button>
          <a v-if="saved" class="bar-btn" :href="`/api/docking/${sid}/zip`">{{ t('web.docking.dl_zip') }}</a>
        </div>

        <h4 class="sub-sec">{{ t('web.docking.sec_3d') }}</h4>
        <div class="set-row">
          <label class="set-label">{{ t('web.docking.pocket_radius') }}</label>
          <Stepper v-model="cutoff" :min="3" :max="6" :step="0.5" />
          <label class="chk"><input v-model="focus" type="checkbox" @change="scheduleReload" /> {{ t('web.docking.focus_chk') }}</label>
          <button class="bar-btn" :disabled="view3dBusy || !scores.length" @click="load3d">{{ t('web.docking.refresh_view') }}</button>
          <button class="bar-btn" :disabled="save3dBusy || !saved" @click="save3d">
            {{ save3dBusy ? t('web.docking.saving_3d') : t('web.docking.save_3d_btn') }}
          </button>
        </div>
        <p class="hint-line">{{ t('web.docking.legend') }}</p>
        <iframe v-if="view3dSrc" :key="view3dKey" class="fig-frame" :src="view3dSrc" title="3d preview" />
        <p v-else class="hint-line">{{ t('web.docking.load_3d_hint') }}</p>
        <p v-if="pocketNames.length" class="hint-line">
          {{ t('web.docking.pocket_list', { n: pocketNames.length, cutoff, list: pocketNames.slice(0, 40).join(', ') }) }}
        </p>
        <template v-if="interRows.length">
          <button class="bar-btn" @click="interOpen = !interOpen">
            {{ interOpen ? t('web.docking.inter_collapse') : t('web.docking.inter_btn', { n: interRows.length }) }}
          </button>
          <div v-if="interOpen" class="tbl-wrap">
            <table class="tbl">
              <thead><tr><th>{{ t('web.docking.th_kind') }}</th><th>{{ t('web.docking.th_residue') }}</th><th>{{ t('web.docking.th_res_atom') }}</th><th>{{ t('web.docking.th_lig_atom') }}</th><th>{{ t('web.docking.th_dist') }}</th></tr></thead>
              <tbody>
                <tr v-for="(r, i) in interRows" :key="i">
                  <td>{{ r.kind === 'hbond' ? t('web.docking.kind_hbond') : r.kind === 'salt' ? t('web.docking.kind_salt') : t('web.docking.kind_hydrophobic') }}</td>
                  <td>{{ (r.chain ? r.chain + '/' : '') + String(r.resn) + String(r.resi) }}</td>
                  <td>{{ r.res_atom }}</td><td>{{ r.lig_atom }}</td><td>{{ fmt(r.dist) }}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </template>
        <p v-else class="hint-line">{{ t('web.docking.no_interactions') }}</p>
      </template>
      <p v-else class="hint-line">{{ t('web.docking.archive_idle') }}</p>
    </section>
  </div>
</template>

<style scoped>
/* 自带受体上传:隐藏原生 file input(原生控件不入整体风格),整行 label 可点击 */
.upload-btn input[type="file"] {
  display: none;
}
.upload-row {
  display: flex;
  align-items: baseline;
  gap: 10px;
  flex-wrap: wrap;
}
</style>

