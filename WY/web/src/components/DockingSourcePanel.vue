<script setup lang="ts">
/** Step 1 数据源与候选基因(Phase 4 D4 自 DockingView 拆出):PPI/富集缓存/
 * 手动基因三来源切换、候选指标与结构筛选。s1 与结构筛选结果由壳层持有
 * (genes 计算、Step 2 药名语境与 Step 6 归档前后缀都要读),本面板只做
 * 采集与展示,经三个事件回传:
 * - picked:来源就绪(标准选择 autoLigand=true / 手动基因提交 false),壳层
 *   写 s1、清下游并按需自动获取配体(Step 2 状态在壳层);
 * - cleared:三态切换,壳层清当前来源与下游;
 * - checked:结构筛选完成(job onDone 或不做筛选的直接路径)。
 * expose 供壳层桥接:manualSession(Step 2 无药名语境)、ligCtx(配体代际
 * 守卫 P0-B 的来源上下文签名)、busy(busyLock 聚合)、allowAf(受体候选
 * 请求要用,勾选框在本面板)。
 * 配体自动获取原本内联在 selectSource 里跨组调用 fetchLigand,拆分后改 emit。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type DockingSourceSummary, type DockingStructCheck } from '../api'
import { useJob } from '../composables/useJob'

const props = defineProps<{
  /** 壳层聚合的在途锁(含本面板 busy):锁住来源下拉,防在途结果落错上下文 */
  busyLock: boolean
  s1: DockingSourceSummary | null
  structResult: DockingStructCheck | null
  ensureSession: () => Promise<string>
  flash: (msg: string) => void
}>()
const emit = defineEmits<{
  (e: 'picked', payload: DockingSourceSummary & { autoLigand: boolean }): void
  (e: 'cleared'): void
  (e: 'checked', r: DockingStructCheck): void
}>()
const { t } = useI18n()

const srcMode = ref<'ppi' | 'enr'>('ppi')
const ppiDirs = ref<{ name: string; display: string }[]>([])
const enrFiles = ref<{ name: string; display: string }[]>([])
const selSource = ref('')
const manualGenes = ref('')
// 手动基因对接专属入口:manual_* 桶不再混入 PPI/富集下拉,由此处接管
const manualMode = ref(false)
const manualDirs = ref<{ name: string; n_genes: number; mtime: string; genes: string[] }[]>([])
const manualPick = ref('')
const useStruct = ref(true)
const allowAf = ref(false)

// 手动会话:Step 2 隐藏"按药名查 PubChem"(无药名语境;CID 直查与 SMILES 保留)
const isManualSession = computed(() => manualMode.value && !props.s1?.drug)

const structJob = useJob((result) => {
  emit('checked', result as unknown as DockingStructCheck)
})

/** 完成徽标与候选计数:与壳层 genes 计算同构(structResult 优先回退 s1 候选) */
const genesCount = computed(() =>
  props.structResult?.genes.length ?? props.s1?.candidates.length ?? 0)

async function loadSources(): Promise<void> {
  const s = await api.dockingSources()
  ppiDirs.value = s.ppi_dirs
  enrFiles.value = s.enr_files
  manualDirs.value = s.manual_dirs ?? []
}

async function selectSource(name: string): Promise<void> {
  const sid = await props.ensureSession()
  const r = await api.dockingSource(sid, { mode: srcMode.value, name })
  // 下游清场与配体自动获取由壳层接手(Step 2 状态在壳层)
  emit('picked', { drug: r.drug, disease: r.disease, hub: r.hub_count, enr: r.enr_count,
                   candidates: r.candidates, autoLigand: true })
}

function pickManualFile(name: string): void {
  const d = manualDirs.value.find((x) => x.name === name)
  if (d?.genes?.length) {
    manualGenes.value = d.genes.join(', ')
    props.flash(t('web.docking.manual_loaded', { n: d.genes.length }))
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
    emit('cleared')
    if (!manualDirs.value.length) void loadSources()
    return
  }
  manualMode.value = false
  if (srcMode.value === mode) return
  srcMode.value = mode
  emit('cleared')
}

async function useManual(): Promise<void> {
  if (!manualGenes.value.trim()) return
  const sid = await props.ensureSession()
  const r = await api.dockingSource(sid, { mode: 'manual', manual_genes: manualGenes.value })
  // 手动路径不动已填的配体查询串(与原 useManual 一致,autoLigand=false)
  emit('picked', { drug: r.drug, disease: r.disease, hub: r.hub_count, enr: r.enr_count,
                   candidates: r.candidates, autoLigand: false })
}

async function checkStructures(): Promise<void> {
  const sid = await props.ensureSession()
  const body = { use_struct: useStruct.value, allow_af: allowAf.value }
  if (!useStruct.value) {
    // 不做结构筛选:候选直接作为最终靶点,无 job
    const r = await api.dockingCheckStructures(sid, body)
    emit('checked', (r.direct ?? {}) as unknown as DockingStructCheck)
    return
  }
  await structJob.run(api.dockingCheckStructures(sid, body) as Promise<{ job_id: string }>,
    { doneText: t('web.docking.check_done_text') })
}

const candidatePreview = computed(() => {
  const list = props.s1?.candidates ?? []
  return list.slice(0, 60).join(', ') + (list.length > 60 ? ' …' : '')
})

onMounted(async () => {
  try {
    await loadSources()
  } catch (e) {
    props.flash(e instanceof Error ? e.message : String(e))
  }
})

/* 壳层桥接位,见文件头 expose 说明 */
defineExpose({
  get manualSession(): boolean { return isManualSession.value },
  get ligCtx(): string {
    return `${srcMode.value}|${selSource.value}|${manualMode.value ? 'manual' : 'auto'}`
  },
  get busy(): boolean { return structJob.running.value },
  get allowAf(): boolean { return allowAf.value }
})
</script>

<template>
  <section class="step" :class="{ done: genesCount > 0 }">
    <h3 class="step-title">{{ t('web.docking.s1_title') }} <span v-if="genesCount > 0" class="ok-mark">{{ t('web.docking.s1_genes_ok', { n: genesCount }) }}</span></h3>

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
</template>
