<script setup lang="ts">
/** ML counts 模式分组选择器:与矩阵模式同外壳、不同发动机(原因见
 * core.expression.counts_picker_fields):正则匹配对象是 counts 列名 →
 * 建议物 = 取值对应列组的最长公共前缀。填入 emit('apply', case, ctrl)
 * 交壳层写 casePat/ctrlPat(值为转义后前缀)。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import type { MlCountsField } from '../api'

const props = defineProps<{
  fields: MlCountsField[]
  /** counts 列名全集:算 spill 用(建议前缀波及组外列,如 D5_20 与 D5_20B 并存) */
  allCols: string[]
  /** 数据上下文 GSE:切换即清本选择器(拆分前由壳层 watch(gse) 统一清) */
  gse: string
  unusable: boolean
  flash: (msg: string) => void
}>()
const emit = defineEmits<{ apply: [case_: string, ctrl: string] }>()
const { t } = useI18n()

const field = ref('')
const caseVal = ref('')
const ctrlVal = ref('')
const current = computed(() => props.fields.find((f) => f.name === field.value) ?? null)
const caseGroup = computed(() =>
  current.value?.groups.find((g) => g.value === caseVal.value) ?? null)
const ctrlGroup = computed(() =>
  current.value?.groups.find((g) => g.value === ctrlVal.value) ?? null)
const spill = computed(() => {
  // 建议前缀会波及组外列(如 D5_20 与 D5_20B 并存)→ 标红,建议手调
  const g = caseGroup.value
  if (!g || !g.suggest) return 0
  return props.allCols.filter((c) => c.includes(g.suggest)).length - g.columns.length
})
const valid = computed(() => {
  if (!caseVal.value || !ctrlVal.value) return null
  if (caseVal.value === ctrlVal.value) return { ok: false, reason: 'same' }
  // 组内列名无公共前缀(如 run accession 混名)→ 建议为空,填入会得到空正则:
  // 空正则匹配一切,canRun 也会因正则为空而锁死,必须在此拦下(QA 2026-11-09)
  if (!caseGroup.value?.suggest || !ctrlGroup.value?.suggest) return { ok: false, reason: 'no_lcp' }
  if ((caseGroup.value?.columns.length ?? 0) < 3
      || (ctrlGroup.value?.columns.length ?? 0) < 3) return { ok: false, reason: 'small' }
  if (spill.value > 0) return { ok: false, reason: 'spill' }
  return { ok: true, reason: '' }
})

watch(() => props.gse, () => {
  field.value = ''
  caseVal.value = ''
  ctrlVal.value = ''
})

function pickField(name: string): void {
  field.value = name
  const f = props.fields.find((x) => x.name === name)
  caseVal.value = f?.suggest_case ?? ''
  ctrlVal.value = f?.suggest_control ?? ''
}

/* 填入成功的按钮旁短反馈:2s 后消失(flash 在页顶,填按钮在页腰,容易两头都看不到) */
const filled = ref(false)
let filledTimer = 0

function apply(): void {
  if (!valid.value?.ok) return
  const esc = (v: string) => v.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')
  emit('apply', esc(caseGroup.value?.suggest ?? ''), esc(ctrlGroup.value?.suggest ?? ''))
  props.flash(t('web.ml.gp_filled'))
  filled.value = true
  window.clearTimeout(filledTimer)
  filledTimer = window.setTimeout(() => { filled.value = false }, 2000)
}
</script>

<template>
  <div class="set-row">
    <label class="set-label">{{ t('web.ml.gp_field_label') }}</label>
    <select class="sel inline" v-model="field" :disabled="unusable"
            @change="pickField(field)">
      <option value="" disabled>{{ t('web.ml.gp_field_ph') }}</option>
      <option v-for="f in fields" :key="f.name" :value="f.name">
        {{ f.label }} ({{ f.groups.length }} 组)
      </option>
    </select>
  </div>
  <template v-if="current">
    <div class="set-row">
      <label class="set-label">{{ t('web.ml.gp_case') }}</label>
      <select class="sel inline" v-model="caseVal" :disabled="unusable">
        <option value="" disabled>{{ t('web.ml.gp_value_ph') }}</option>
        <option v-for="g in current.groups" :key="g.value" :value="g.value">
          {{ g.value.slice(0, 60) }} → {{ g.columns.length }} 列 · {{ g.suggest ? `建议 ${g.suggest}` : t('web.ml.gp_no_lcp') }}
        </option>
      </select>
    </div>
    <div class="set-row">
      <label class="set-label">{{ t('web.ml.gp_ctrl') }}</label>
      <select class="sel inline" v-model="ctrlVal" :disabled="unusable">
        <option value="" disabled>{{ t('web.ml.gp_value_ph') }}</option>
        <option v-for="g in current.groups" :key="g.value" :value="g.value">
          {{ g.value.slice(0, 60) }} → {{ g.columns.length }} 列 · {{ g.suggest ? `建议 ${g.suggest}` : t('web.ml.gp_no_lcp') }}
        </option>
      </select>
    </div>
    <p v-if="caseVal && ctrlVal" class="hint-line"
       :class="{ good: valid?.ok, stale: !valid?.ok }">
      {{ valid?.ok
        ? t('web.ml.gp_covers', { n: caseGroup?.columns.length, m: ctrlGroup?.columns.length })
        : t('web.ml.gp_validate_fail',
            { reason: t('web.ml.gp_reason_' + valid?.reason) }) }}
    </p>
    <button class="bar-btn" :disabled="!valid?.ok || unusable" @click="apply">
      {{ t('web.ml.gp_fill') }}
    </button>
    <span v-if="filled" class="state on">{{ t('web.ml.gp_filled') }}</span>
  </template>
</template>
