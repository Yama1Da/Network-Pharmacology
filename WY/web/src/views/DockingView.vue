<script setup lang="ts">
/** 分子对接(阶段 3b,Molecular_Docking 的 TS 迁移):六步向导。
 * 会话状态存服务端(core.docking_session),前端持 sid;网络重步骤走 job 轮询;
 * 3D 预览为服务端 py3Dmol 渲染的 HTML,iframe 内嵌。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type Row } from '../api'
import Stepper from '../components/Stepper.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'

const { t } = useI18n()

/* ---------------- 会话 ---------------- */

const sid = ref('')
const dockError = ref('')

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
const s1 = ref<{ drug: string; disease: string; hub: number; enr: number; candidates: string[] } | null>(null)
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

async function useManual(): Promise<void> {
  sid.value = await ensureSession()
  const r = await api.dockingSource(sid.value, { mode: 'manual', manual_genes: manualGenes.value })
  s1.value = { drug: r.drug, disease: r.disease, hub: r.hub_count, enr: r.enr_count, candidates: r.candidates }
  structResult.value = null
  resetDownstream()
}

/** 数据源变更后清空前端本地的下游视图状态(服务端会话已由端点同步失效)。 */
function resetDownstream(): void {
  lig.value = null
  receptor.value = null
  candRows.value = []
  candError.value = ''
  pdbSel.value = ''
  gene.value = ''
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

async function fetchLigand(): Promise<void> {
  sid.value = await ensureSession()
  ligBusy.value = true
  try {
    lig.value = ligSource.value === 'pubchem'
      ? await api.dockingLigand(sid.value, { source: 'pubchem', query: ligQuery.value })
      : await api.dockingLigand(sid.value, { source: 'smiles', smiles: ligSmiles.value })
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
const recJob = useJob((result) => {
  receptor.value = result as unknown as { pdb_id: string; receptor_source: string }
  flash(t('web.docking.receptor_ready'))
})

async function loadCandidates(): Promise<void> {
  if (!gene.value) return
  sid.value = await ensureSession()
  pdbSel.value = ''
  receptor.value = null
  box.value = null
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
  await recJob.run(api.dockingReceptorPrepare(sid.value, { gene: gene.value, pdb_id: pdbSel.value }),
    { doneText: t('web.docking.receptor_ready') })
}

const selFlag = computed(() => candRows.value.find((r) => r.pdb_id === pdbSel.value)?.flag ?? '')

/* ---------------- Step 4 对接盒 ---------------- */

const box = ref<{ center: number[]; size: number[] } | null>(null)
const boxMeta = ref<{ clamped: boolean; inferred_ligands: string[] | null; residues: string[] } | null>(null)
const boxBusy = ref(false)

async function inferBox(): Promise<void> {
  sid.value = await ensureSession()
  boxBusy.value = true
  try {
    const r = await api.dockingBox(sid.value, { infer: true })
    box.value = r.box
    boxMeta.value = { clamped: r.clamped, inferred_ligands: r.inferred_ligands, residues: r.residues }
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
  } catch (e) {
    flash(e instanceof Error ? e.message : String(e))
  } finally {
    boxBusy.value = false
  }
}

function editCenter(i: number, e: Event): void {
  const v = Number((e.target as HTMLInputElement).value)
  if (!box.value) return
  box.value = { center: box.value.center.map((c, j) => (j === i ? v : c)),
                size: box.value.size }
}
function editSize(i: number, e: Event): void {
  const v = Number((e.target as HTMLInputElement).value)
  if (!box.value) return
  box.value = { center: box.value.center,
                size: box.value.size.map((s, j) => (j === i ? v : s)) }
}

const fmt = (v: unknown, digits = 3): string =>
  typeof v === 'number' ? v.toFixed(digits) : (v == null ? '—' : String(v))

/* ---------------- Step 5 运行对接 ---------------- */

const runJob = useJob((result) => {
  const r = result as { log_tail: string; scores: Row[] }
  dockState.value = { result: { log_tail: r.log_tail, scores: r.scores } }
})
const dockState = ref<{ result: { log_tail: string; scores: Row[] } } | null>(null)

const canRun = computed(() => !!(
  lig.value?.ok && receptor.value && box.value))

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
    const r = await api.dockingSave3d(sid.value, { pose: poseSel.value, caption })
    flash(t('web.docking.saved_3d', { rel: r.artifact_rel,
                                       offline: r.offline_ok ? t('web.docking.offline_ok') : t('web.docking.offline_no') }))
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
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <!-- ============ Step 1 ============ -->
    <section class="step" :class="{ done: stepDone.s1 }">
      <h3 class="step-title">{{ t('web.docking.s1_title') }} <span v-if="stepDone.s1" class="ok-mark">{{ t('web.docking.s1_genes_ok', { n: genes.length }) }}</span></h3>
      <div class="lv-group">
        <button class="lv" :class="{ on: srcMode === 'ppi' }" @click="srcMode = 'ppi'">{{ t('web.docking.src_ppi') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: srcMode === 'enr' }" @click="srcMode = 'enr'">{{ t('web.docking.src_enr') }}</button>
      </div>
      <p v-if="srcMode === 'ppi' && !ppiDirs.length" class="hint-line">{{ t('web.docking.no_ppi_dirs') }}</p>
      <select v-if="srcMode === 'ppi' && ppiDirs.length" class="sel" v-model="selSource" @change="selectSource(selSource)">
        <option value="" disabled>{{ t('web.docking.sel_dir_ph') }}</option>
        <option v-for="d in ppiDirs" :key="d.name" :value="d.name">{{ d.display }}</option>
      </select>
      <select v-if="srcMode === 'enr' && enrFiles.length" class="sel" v-model="selSource" @change="selectSource(selSource)">
        <option value="" disabled>{{ t('web.docking.sel_enr_ph') }}</option>
        <option v-for="f in enrFiles" :key="f.name" :value="f.name">{{ f.display }}</option>
      </select>

      <template v-if="s1">
        <div class="metrics three">
          <div class="metric"><b>{{ s1.hub }}</b><span>{{ t('web.docking.m_hub') }}</span></div>
          <div class="metric"><b>{{ s1.enr }}</b><span>{{ t('web.docking.m_enr') }}</span></div>
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

      <div v-if="!s1 && !ppiDirs.length && !enrFiles.length" class="manual">
        <p class="hint-line">{{ t('web.docking.manual_hint') }}</p>
        <textarea v-model="manualGenes" rows="2" class="ta" placeholder="TP53, EGFR, …"></textarea>
        <button class="bar-btn" @click="useManual">{{ t('web.docking.manual_btn') }}</button>
      </div>
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
          <label class="set-label">{{ t('web.docking.lig_label') }}</label>
          <div class="field"><input v-model="ligQuery" :placeholder="s1?.drug || t('web.docking.drug_name_ph')" spellcheck="false" /><i class="line"></i></div>
          <button class="bar-btn" :disabled="ligBusy" @click="fetchLigand">
            {{ ligBusy ? t('web.docking.fetching') : t('web.docking.fetch_ligand') }}
          </button>
        </div>
        <p v-if="!s1?.drug && !lig?.ok" class="hint-line">{{ t('web.docking.no_drug_hint') }}</p>
      </template>
      <template v-else>
        <textarea v-model="ligSmiles" rows="2" class="ta" :placeholder="t('web.docking.smiles_ph')"></textarea>
        <button class="bar-btn" :disabled="ligBusy || !ligSmiles.trim()" @click="fetchLigand">
          {{ ligBusy ? t('web.docking.generating') : t('web.docking.gen_conformer') }}
        </button>
      </template>
      <p v-if="lig?.ok" class="hint-line good">{{ t('web.docking.lig_ok') }}</p>
      <p v-if="lig && !lig.ok" class="stale">{{ lig.error }}</p>
    </section>

    <!-- ============ Step 3 ============ -->
    <section class="step" :class="{ done: stepDone.s3 }">
      <h3 class="step-title">{{ t('web.docking.s3_title') }} <span v-if="stepDone.s3" class="ok-mark">✓ {{ receptor?.pdb_id }}</span></h3>
      <template v-if="genes.length">
        <select class="sel" v-model="gene" @change="loadCandidates">
          <option value="" disabled>{{ t('web.docking.sel_gene_ph') }}</option>
          <option v-for="g in genes" :key="g" :value="g">{{ g }}</option>
        </select>
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
          <div v-if="recJob.running.value" class="progress">
            <div class="progress-track">
              <div class="progress-fill" :style="{ width: `${((recJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
            </div>
            <p class="progress-text">{{ recJob.text.value }}</p>
          </div>
          <p v-if="receptor?.receptor_source === 'alphafold'" class="stale">
            {{ t('web.docking.af_warn', { pdb: receptor.pdb_id.split(':')[1] }) }}
          </p>
        </template>
        <p v-else-if="gene && !candError && !recJob.running.value" class="hint-line">
          {{ t('web.docking.rcsb_empty') }}
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
      <p v-if="boxMeta?.clamped" class="stale">{{ t('web.docking.clamped_note') }}</p>
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
            <iframe class="fig-frame" :src="`/api/docking/${sid}/box-3d`" title="box preview" />
          </div>
        </div>
      </template>
    </section>

    <!-- ============ Step 5 ============ -->
    <section class="step" :class="{ done: stepDone.s5 }">
      <h3 class="step-title">{{ t('web.docking.s5_title') }} <span v-if="stepDone.s5" class="ok-mark">✓</span></h3>
      <p v-if="canRun" class="hint-line">{{ t('web.docking.run_params') }}</p>
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
          <button class="bar-btn strong" :disabled="saveBusy" @click="saveResult">
            {{ saveBusy ? t('web.docking.archiving') : t('web.docking.archive_btn') }}
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
        <iframe v-if="view3dSrc" class="fig-frame" :src="view3dSrc" title="3d preview" />
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

