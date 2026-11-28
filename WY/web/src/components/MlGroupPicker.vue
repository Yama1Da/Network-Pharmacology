<script setup lang="ts">
/** ML 矩阵模式分组选择器(总览的三下拉:字段 → case/control 取值)。
 * 生成的正则为「列取值全值」的转义 —— 与 assign_groups 整行拼接匹配同源,
 * protocol 列复述队列描述造成的污染从源头不存在;重叠对由后端 overlaps 显式给出。
 * 填入 emit('apply', case, ctrl) 交壳层写 casePat/ctrlPat(值为转义后正则)。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import type { MlGroupField } from '../api'

const props = defineProps<{
  fields: MlGroupField[]
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
const caseN = computed(() =>
  current.value?.values.find((v) => v.value === caseVal.value)?.n ?? null)
const ctrlN = computed(() =>
  current.value?.values.find((v) => v.value === ctrlVal.value)?.n ?? null)
const overlap = computed(() => {
  const f = current.value
  if (!f || !caseVal.value || !ctrlVal.value) return null
  return f.overlaps.find((o) => (o.a === caseVal.value && o.b === ctrlVal.value)
                             || (o.a === ctrlVal.value && o.b === caseVal.value)) ?? null
})
const valid = computed(() => {
  if (!caseVal.value || !ctrlVal.value) return null
  if (caseVal.value === ctrlVal.value) return { ok: false, reason: 'same' }
  if (overlap.value) return { ok: false, reason: 'overlap' }
  if ((caseN.value ?? 0) === 0 || (ctrlN.value ?? 0) === 0) return { ok: false, reason: 'zero' }
  if ((caseN.value ?? 0) < 3 || (ctrlN.value ?? 0) < 3) return { ok: false, reason: 'small' }
  return { ok: true, reason: '' }
})
/* 疑似纵向字段(取值含治疗前后/周次等):病例组含配对样本时提示用户收紧到时间点 */
const pairedField = computed(() => {
  for (const f of props.fields) {
    if (f.name === field.value) continue
    if (f.values.some((v) =>
      /before|after|week|pre-?treatment|post-?treatment|day\s*\d/i.test(v.value)))
      return f.label
  }
  return ''
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
  emit('apply', esc(caseVal.value), esc(ctrlVal.value))
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
        {{ f.label }} ({{ f.values.length }})
      </option>
    </select>
  </div>
  <template v-if="current">
    <div class="set-row">
      <label class="set-label">{{ t('web.ml.gp_case') }}</label>
      <select class="sel inline" v-model="caseVal" :disabled="unusable">
        <option value="" disabled>{{ t('web.ml.gp_value_ph') }}</option>
        <option v-for="v in current.values" :key="v.value" :value="v.value">
          {{ v.value }} ({{ v.n }})
        </option>
      </select>
    </div>
    <div class="set-row">
      <label class="set-label">{{ t('web.ml.gp_ctrl') }}</label>
      <select class="sel inline" v-model="ctrlVal" :disabled="unusable">
        <option value="" disabled>{{ t('web.ml.gp_value_ph') }}</option>
        <option v-for="v in current.values" :key="v.value" :value="v.value">
          {{ v.value }} ({{ v.n }})
        </option>
      </select>
    </div>
    <p v-if="pairedField" class="stale">
      {{ t('web.ml.gp_paired', { field: pairedField }) }}
    </p>
    <p v-if="caseVal && ctrlVal" class="hint-line"
       :class="{ good: valid?.ok, stale: !valid?.ok }">
      {{ valid?.ok
        ? t('web.ml.gp_validate_ok', { n: caseN, m: ctrlN })
        : t('web.ml.gp_validate_fail',
            { reason: t('web.ml.gp_reason_' + valid?.reason) }) }}
    </p>
    <button class="bar-btn" :disabled="!valid?.ok || unusable" @click="apply">
      {{ t('web.ml.gp_fill') }}
    </button>
    <span v-if="filled" class="state on">{{ t('web.ml.gp_filled') }}</span>
  </template>
</template>
