<script setup lang="ts">
/** Step 6 结果与归档(Phase 4 D1 自 DockingView 拆出):构象选择、结论、归档
 * 防重复签名、3D 视图(py3Dmol iframe)与互作表。poseSel 由壳层持有经 v-model
 * 双向——Step 5 得分表行点击与本面板构象下拉共写一个值;受体 pdb_id 与对接
 * 轮次经 props 参与归档签名,任一变化重新允许归档。expose clear() 供壳层
 * resetSteps4to6 / resetDownstream 显式编排清场。 */
import { computed, onBeforeUnmount, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl, type Row } from '../api'
import Stepper from './Stepper.vue'
import { fmt } from '../utils/format'

const props = defineProps<{
  sid: string
  scores: Row[]
  /** 受体 pdb_id(归档签名成分;受体本身状态在壳层,经受体面板 v-model) */
  pdbId: string
  gene: string
  drug: string
  disease: string
  /** 对接轮次(壳层 runSeq):新一轮结果未归档时签名不同,重新允许归档 */
  runSeq: number
  poseSel: number
  ensureSession: () => Promise<string>
  flash: (msg: string) => void
}>()
const emit = defineEmits<{ (e: 'update:poseSel', v: number): void }>()
const { t } = useI18n()

/* poseSel 的 v-model 桥:构象下拉写回壳层,再由壳层带给 Step 5 得分表高亮 */
const poseModel = computed({
  get: () => props.poseSel,
  set: (v: number) => emit('update:poseSel', v)
})

const conclusion = ref('')
const drugName = ref('')
const diseaseName = ref('')
const saved = ref<{ result_dir: string } | null>(null)
const starred = ref(false)
const saveBusy = ref(false)

/* 归档防重复:签名 = 对接轮次 + 受体 + 构象 + 结论 + 前后缀;归档后按钮置灰,
 * 上述任一内容变化才重新允许归档(想留多版本就改一下结论再归)。 */
const lastSavedSig = ref('')
const archiveSig = computed(() => JSON.stringify(
  [props.runSeq, props.pdbId, props.poseSel,
   conclusion.value, drugName.value, diseaseName.value]))
const archiveDirty = computed(() => archiveSig.value !== lastSavedSig.value)

const cutoff = ref(4.0)
const focus = ref(false)
const view3dKey = ref(0)
const view3dSrc = ref('')
const view3dBusy = ref(false)
const pocketRows = ref<{ resn: string; resi: string; chain: string; dist: number }[]>([])
const pocketNames = computed(() => pocketRows.value.map((r) => `${r.resn}${r.resi}`))
const interRows = ref<Row[]>([])
const interOpen = ref(false)

async function load3d(): Promise<void> {
  if (!props.scores.length || !props.sid) return
  view3dBusy.value = true
  try {
    view3dKey.value++  // 重建 iframe:同参数再点"刷新视图"也能真正重载,而不是静默无反应
    view3dSrc.value = api.dockingView3dUrl(props.sid, props.poseSel, cutoff.value, focus.value)
    const inter = await api.dockingInteractions(props.sid, props.poseSel)
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
onBeforeUnmount(() => window.clearTimeout(reloadTimer))

async function saveResult(): Promise<void> {
  const sid = await props.ensureSession()
  saveBusy.value = true
  try {
    const r = await api.dockingSave(sid, {
      pose_mode: props.poseSel, conclusion: conclusion.value,
      drug: drugName.value, disease: diseaseName.value
    })
    saved.value = { result_dir: r.result_dir }
    lastSavedSig.value = archiveSig.value
    starred.value = false
    props.flash(t('web.docking.archived'))
  } catch (e) {
    props.flash(e instanceof Error ? e.message : String(e))
  } finally {
    saveBusy.value = false
  }
}

async function toggleStar(): Promise<void> {
  if (!saved.value) return
  try {
    starred.value = (await api.dockingStar(props.sid, !starred.value)).starred
  } catch {
    props.flash(t('web.docking.star_fail'))
  }
}

const save3dBusy = ref(false)
async function save3d(): Promise<void> {
  if (!saved.value) return
  save3dBusy.value = true
  try {
    const caption = `NPAP · ${drugName.value || '-'} × ${diseaseName.value || '-'} · ${props.gene} · pose ${props.poseSel}`
    const r = await api.dockingSave3d(props.sid,
      { pose: props.poseSel, caption, cutoff: cutoff.value })
    props.flash(t('web.docking.saved_3d', { rel: r.artifact_rel,
                                             offline: r.offline_ok ? t('web.docking.offline_ok') : t('web.docking.offline_no') }))
    // 保存进结果目录的同时把 HTML 直接下载到本地(对齐旧版"下载交互式 3D 视图")
    const a = document.createElement('a')
    a.href = artifactsUrl(r.artifact_rel, true)
    a.click()
  } catch (e) {
    props.flash(e instanceof Error ? e.message : String(e))
  } finally {
    save3dBusy.value = false
  }
}

/* 壳层编排位:受体/数据源变更后清掉本步旧结果(归档、结论、3D 与互作都绑定
 * 具体受体与具体一轮对接)。drugName/diseaseName 与 lastSavedSig 不清,与原
 * DockingView 的两处 reset 保持一致。 */
defineExpose({
  clear(): void {
    conclusion.value = ''
    saved.value = null
    starred.value = false
    view3dSrc.value = ''
    view3dKey.value++
    interRows.value = []
    pocketRows.value = []
    interOpen.value = false
  }
})
</script>

<template>
  <section class="step" :class="{ done: !!saved }">
    <h3 class="step-title">{{ t('web.docking.s6_title') }}
      <span v-if="saved" class="ok-mark">✓ {{ saved.result_dir.split('/').pop() || saved.result_dir }}</span>
    </h3>
    <template v-if="scores.length">
      <div class="set-row">
        <label class="set-label">{{ t('web.docking.pose_label') }}</label>
        <select class="sel inline" v-model.number="poseModel" @change="scheduleReload">
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
        <div class="field"><input v-model="drugName" :placeholder="drug || t('web.docking.drug_name_ph')" spellcheck="false" /><i class="line"></i></div>
        <div class="field"><input v-model="diseaseName" :placeholder="disease || t('web.docking.disease_name_ph')" spellcheck="false" /><i class="line"></i></div>
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
</template>
