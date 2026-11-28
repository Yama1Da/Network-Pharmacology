<script setup lang="ts">
/** Step 4 对接盒(Phase 4 D2 自 DockingView 拆出):自动推断/手动微调对接盒
 * 坐标与盒预览 iframe。box 由壳层持有经 v-model 双向——壳层 canRun 与
 * Step 4 完成徽标都要读它;boxMeta/boxBusy/boxKey 为组件私有。expose clear()
 * 供壳层 resetSteps4to6 / resetDownstream 显式编排清场。 */
import { ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from '../api'

const props = defineProps<{
  sid: string
  receptor: { pdb_id: string; receptor_source: string } | null
  box: { center: number[]; size: number[] } | null
  ensureSession: () => Promise<string>
  flash: (msg: string) => void
}>()
const emit = defineEmits<{
  (e: 'update:box', v: { center: number[]; size: number[] } | null): void
}>()
const { t } = useI18n()

const boxMeta = ref<{ clamped: boolean; inferred_ligands: string[] | null; residues: string[] } | null>(null)
const boxBusy = ref(false)
// 盒子预览重建计数:iframe 的 src 固定指向 box-3d,微调/重推不重建就永远不刷新
// (QA 2026-11-09 反馈)。递增同时改 src(?v=)与 :key,防缓存且强制重挂载。
const boxKey = ref(0)

async function inferBox(): Promise<void> {
  const sid = await props.ensureSession()
  boxBusy.value = true
  try {
    const r = await api.dockingBox(sid, { infer: true })
    emit('update:box', r.box)
    boxMeta.value = { clamped: r.clamped, inferred_ligands: r.inferred_ligands, residues: r.residues }
    boxKey.value++
  } catch (e) {
    props.flash(e instanceof Error ? e.message : String(e))
  } finally {
    boxBusy.value = false
  }
}

async function applyManualBox(): Promise<void> {
  const sid = await props.ensureSession()
  boxBusy.value = true
  try {
    const r = await api.dockingBox(sid, { infer: false, center: props.box!.center, size: props.box!.size })
    boxMeta.value = { clamped: false, inferred_ligands: null, residues: r.residues }
    boxKey.value++
  } catch (e) {
    props.flash(e instanceof Error ? e.message : String(e))
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
  if (v === null || !props.box) return
  emit('update:box', { center: props.box.center.map((c, j) => (j === i ? v : c)),
                       size: props.box.size })
}
function editSize(i: number, e: Event): void {
  const v = parseBoxInput(e)
  if (v === null || !props.box) return
  emit('update:box', { center: props.box.center,
                       size: props.box.size.map((s, j) => (j === i ? v : s)) })
}

/* 壳层编排位:换受体/换数据源后旧盒子整体作废(盒坐标绑定具体受体)。 */
defineExpose({
  clear(): void {
    emit('update:box', null)
    boxMeta.value = null
  }
})
</script>

<template>
  <section class="step" :class="{ done: !!box }">
    <h3 class="step-title">{{ t('web.docking.s4_title') }} <span v-if="box" class="ok-mark">✓</span></h3>
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
</template>
