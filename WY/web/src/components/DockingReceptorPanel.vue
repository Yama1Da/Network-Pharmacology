<script setup lang="ts">
/** Step 3 受体(Phase 4 D3 自 DockingView 拆出):基因下拉 + 双入路( RCSB 检索
 * 制备 / 自带结构上传),候选表、制备 job 与上传体检报告随组件。receptor 由
 * 壳层持有经 v-model 双向——壳层 canRun、Step 3 完成徽标与对接盒面板都要读。
 * 原本四处内联的 resetSteps4to6() 改为 emit('reset-downstream') 交壳层编排
 * (清 Step 4-6 的状态如今分属对接盒/归档面板组件)。expose:clear() 供壳层
 * resetDownstream 编排;busy 供壳层聚合 busyLock;gene 供壳层喂归档面板
 * (结论占位与 3D 图注要用)。 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from '../api'
import { useJob } from '../composables/useJob'

const props = defineProps<{
  /** 候选基因(壳层 genes:结构筛选结果 ?? 来源候选) */
  genes: string[]
  /** 允许 AlphaFold 模型(勾选框在 Step 1,值经壳层桥接) */
  allowAf: boolean
  receptor: { pdb_id: string; receptor_source: string } | null
  ensureSession: () => Promise<string>
  flash: (msg: string) => void
}>()
const emit = defineEmits<{
  (e: 'update:receptor', v: { pdb_id: string; receptor_source: string } | null): void
  /** 本步变化令 Step 4-6 的旧结果全部失效(换基因/换受体/换入路/换文件) */
  (e: 'reset-downstream'): void
}>()
const { t } = useI18n()

const gene = ref('')
const candRows = ref<{ pdb_id: string; resolution: string; method: string; title: string; source: string; flag: string; note: string }[]>([])
const candError = ref('')
const pdbSel = ref('')
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
  emit('update:receptor', r)
  uploadReport.value = r.report ?? null
  emit('reset-downstream')
  props.flash(t(r.report ? 'web.docking.upload_ready' : 'web.docking.receptor_ready'))
})

const recFailLabel = computed(() =>
  recSource.value === 'upload' ? t('web.docking.upload_fail') : t('web.docking.prep_fail'))

/** 切换受体来源:清 receptor 与下游(两条入路的结果互斥),候选表保留便于切回 */
function setRecSource(mode: 'rcsb' | 'upload'): void {
  if (recSource.value === mode || recJob.running.value) return
  recSource.value = mode
  emit('update:receptor', null)
  pdbSel.value = ''
  pickedFile.value = null
  uploadReport.value = null
  emit('reset-downstream')
}

function onFilePicked(e: Event): void {
  const input = e.target as HTMLInputElement
  const f = input.files?.[0] ?? null
  input.value = '' // 允许重复选择同一个文件重新上传
  if (!f) return
  pickedFile.value = f
  uploadReport.value = null
  if (props.receptor) {  // 换文件即换结构:旧受体作废
    emit('update:receptor', null)
    emit('reset-downstream')
  }
}

async function uploadReceptor(): Promise<void> {
  if (!gene.value || !pickedFile.value) return
  const sid = await props.ensureSession()
  // 代际守卫(P0-B):基因切换或切回 RCSB 分支都视为上下文失效
  const g0 = gene.value
  await recJob.run(api.dockingReceptorUpload(sid, gene.value, pickedFile.value),
    { doneText: t('web.docking.upload_ready'),
      stale: () => gene.value !== g0 || recSource.value !== 'upload' })
}

async function loadCandidates(): Promise<void> {
  if (!gene.value) return
  const sid = await props.ensureSession()
  pdbSel.value = ''
  emit('update:receptor', null)
  emit('reset-downstream')
  try {
    const r = await api.dockingReceptorCandidates(sid, { gene: gene.value, allow_af: props.allowAf })
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
  const sid = await props.ensureSession()
  // 代际守卫(P0-B):发起时捕获基因快照,期间切换则受体结果作废,
  // 不落 Step 3 展示(服务端另有写前核对兜底)
  const g0 = gene.value
  await recJob.run(api.dockingReceptorPrepare(sid, { gene: gene.value, pdb_id: pdbSel.value }),
    { doneText: t('web.docking.receptor_ready'),
      stale: () => gene.value !== g0 })
}

const selFlag = computed(() => candRows.value.find((r) => r.pdb_id === pdbSel.value)?.flag ?? '')

/* 壳层编排位:换数据源后清空本步全部状态(receptor 由壳层经 v-model 直清)。
 * recSource 不清——与原 DockingView 的 resetDownstream 保持一致。 */
defineExpose({
  clear(): void {
    candRows.value = []
    candError.value = ''
    pdbSel.value = ''
    gene.value = ''
    pickedFile.value = null
    uploadReport.value = null
  },
  /** 壳层 busyLock 聚合位:制备/上传在途时锁住来源切换入口(P0-B) */
  get busy(): boolean { return recJob.running.value },
  /** 壳层桥接位:归档面板的结论占位与 3D 图注需要当前基因 */
  get gene(): string { return gene.value }
})
</script>

<template>
  <section class="step" :class="{ done: !!receptor }">
    <h3 class="step-title">{{ t('web.docking.s3_title') }} <span v-if="receptor" class="ok-mark">✓ {{ receptor.pdb_id }}</span></h3>
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
        <p v-if="!gene" class="hint-line">{{ t('web.docking.upload_need_gene') }}</p>
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
