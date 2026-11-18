<script setup lang="ts">
/** 协同预测页(二期):药物对表单 + 页级参数(N/置换轮数,Settings 作初始默认,页内覆盖)
 * + jobs 进度(stage 码映射文案,fraction 由核心按零分布轮数推进)+ 契约 synergy.v1.1 结果渲染。
 * 措辞红线(清单 §1.3):geometry 三枚举只述几何事实,禁译"协同/拮抗";|z_sep|<1 附弱信号提示。
 * 深色适配:仅家族类 + 布局 scoped 样式,色值全走 var(--c-*),无硬编码。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SynergyContract, type SynergyProx } from '../api'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import { fmt } from '../utils/format'
import { takePrefill } from '../store'

const { t } = useI18n()
const { notice, flash } = useNotice(3200)

const drugA = ref('')
const drugB = ref('')
const disease = ref('')
const drugOptions = ref<string[]>([])
const diseaseOptions = ref<string[]>([])
const N_OPTIONS: (number | null)[] = [30, 50, 100, null]      // null = 全量
const PERM_OPTIONS = [50, 100, 200, 400]
const topN = ref<number | null>(100)
const nPerm = ref<number | null>(200)
const nPermDefault = ref<number>(200)                          // Settings 初始默认,展示兜底用

const result = ref<SynergyContract | null>(null)

const job = useJob((r) => {
  result.value = r as unknown as SynergyContract
})

onMounted(async () => {
  const pf = takePrefill()
  if (pf) {
    drugA.value = pf.drug
    disease.value = pf.disease
  }
  try {
    const lib = await api.synergyLibrary()
    drugOptions.value = lib.drugs
    diseaseOptions.value = lib.diseases
  } catch {
    /* 本地库清单不可达时输入框照常手输 */
  }
  try {
    const s = await api.settings()
    topN.value = s.synergy_top_n ?? 100
    nPerm.value = s.synergy_n_perm ?? 200
    nPermDefault.value = s.synergy_n_perm ?? 200
  } catch {
    /* 设置不可达时保持页面内建默认 */
  }
})

function run(): void {
  if (!drugA.value.trim() || !drugB.value.trim()) {
    flash(t('web.synergy.err_required'))
    return
  }
  void job.run(api.synergyAnalyze({
    drug_a: drugA.value.trim(),
    drug_b: drugB.value.trim(),
    disease: disease.value.trim(),
    top_n: topN.value,
    n_perm: nPerm.value
  }))
}

/* 后端错误按稳定令牌映射 web.synergy.*(api/synergy.py 约定);未识别令牌原样透出 */
const errorMessage = computed(() => {
  const e = job.error.value
  if (!e) return ''
  const map: [string, string][] = [
    ['too_few_targets_a', t('web.synergy.err_too_few_a')],
    ['too_few_targets_b', t('web.synergy.err_too_few_b')],
    ['no_targets_a', t('web.synergy.err_no_targets_a')],
    ['no_targets_b', t('web.synergy.err_no_targets_b')],
    ['no_targets_disease', t('web.synergy.err_no_targets_disease')],
    ['InteractomeNotBuiltError', t('web.synergy.err_not_built')]
  ]
  for (const [token, text] of map) {
    if (e.includes(token)) return text
  }
  return e
})

/* 进度文案:优先 stage 码映射(核心只发 fraction,不带文案),兜底后端 text */
const stageText = computed(() => {
  const stage = String((job.extra.value as Record<string, unknown>).stage ?? '')
  if (stage === 'targets') return t('web.synergy.stage_targets')
  if (stage === 'separation') return t('web.synergy.stage_separation')
  if (stage === 'proximity') return t('web.synergy.stage_proximity')
  return job.text.value
})

const labelText = computed(() => {
  const m: Record<string, string> = {
    overlap: t('web.synergy.label_overlap'),
    separated: t('web.synergy.label_separated'),
    tie: t('web.synergy.label_tie')
  }
  return result.value ? (m[result.value.geometry_label] ?? result.value.geometry_label) : ''
})

/* 弱信号(2026-11-18 实证:宽靶点集/噪声线附近标签随符号翻转):零分布退化或 |z_sep|<1 */
const weakHint = computed(() => {
  const s = result.value?.separation
  if (!s) return false
  return s.null.z === null || Math.abs(s.null.z) < 1
})

/* proximity 三行渲染:ok=true 显示数值;否则透出后端错误令牌(不伪装成数值) */
const proxRows = computed(() => {
  const p = result.value?.proximity
  if (!p || !result.value?.disease) return []
  const rows: { key: string; label: string; v: SynergyProx | null }[] = [
    { key: 'z_a', label: t('web.synergy.res_prox_z_a'), v: p.z_a },
    { key: 'z_b', label: t('web.synergy.res_prox_z_b'), v: p.z_b },
    { key: 'z_ab', label: t('web.synergy.res_prox_z_ab'), v: p.z_ab }
  ]
  return rows
})

const f = (x: number | null | undefined, digits = 4): string =>
  x === null || x === undefined ? t('web.synergy.value_na') : fmt(x, digits)

const topText = (n: number | null): string =>
  n === null ? t('web.synergy.param_n_full') : String(n)

const permText = (n: number | null): string =>
  n === null ? `${nPermDefault.value}${t('web.synergy.res_default_suffix')}` : String(n)
</script>

<template>
  <p class="kicker s1">{{ t('web.synergy.kicker') }}</p>
  <h1 class="title s2">{{ t('web.synergy.title') }}</h1>
  <p class="sub s3">{{ t('web.synergy.sub') }}</p>

  <div class="panel s4">
    <datalist id="synergy-drug-options">
      <option v-for="d in drugOptions" :key="d" :value="d"></option>
    </datalist>
    <datalist id="synergy-disease-options">
      <option v-for="d in diseaseOptions" :key="d" :value="d"></option>
    </datalist>

    <div class="names-form">
      <label class="field">
        <input v-model="drugA" spellcheck="false" list="synergy-drug-options"
               autocomplete="off" :placeholder="t('web.synergy.form_drug_a')" />
        <i class="line"></i>
      </label>
      <span class="cross">×</span>
      <label class="field">
        <input v-model="drugB" spellcheck="false" list="synergy-drug-options"
               autocomplete="off" :placeholder="t('web.synergy.form_drug_b')" />
        <i class="line"></i>
      </label>
    </div>

    <div class="disease-block">
      <span class="ctl-label">{{ t('web.synergy.form_disease') }}</span>
      <label class="field">
        <input v-model="disease" spellcheck="false" list="synergy-disease-options"
               autocomplete="off" :placeholder="t('web.synergy.ph_disease')" />
        <i class="line"></i>
      </label>
    </div>

    <div class="param-row">
      <span class="ctl-label">{{ t('web.synergy.param_n_label') }}</span>
      <span class="lv-group">
        <button v-for="opt in N_OPTIONS" :key="String(opt)" class="lv" type="button"
                :class="{ on: topN === opt }" @click="topN = opt">
          {{ opt === null ? t('web.synergy.param_n_full') : opt }}
        </button>
      </span>
      <span class="hint">{{ t('web.synergy.param_n_hint') }}</span>
    </div>
    <div class="param-row">
      <span class="ctl-label">{{ t('web.synergy.param_perm_label') }}</span>
      <span class="lv-group">
        <button v-for="opt in PERM_OPTIONS" :key="opt" class="lv" type="button"
                :class="{ on: nPerm === opt }" @click="nPerm = opt">{{ opt }}</button>
      </span>
      <span class="hint">{{ t('web.synergy.param_perm_hint') }}</span>
    </div>

    <div class="xlinks">
      <button class="bar-btn strong" type="button" :disabled="job.running.value" @click="run">
        {{ t('web.synergy.run_btn') }}
      </button>
    </div>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ stageText || t('web.synergy.stage_targets') }}</p>
    </div>
    <p v-if="errorMessage" class="notice err">{{ errorMessage }}</p>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>

  <div v-if="result" class="panel s5">
    <p class="res-kicker">{{ t('web.synergy.res_title') }} · {{ result.schema }}</p>

    <div class="geometry">
      <span class="geo-label">{{ labelText }}</span>
      <span class="geo-basis">{{ result.label_basis }}</span>
      <span v-if="weakHint" class="weak">{{ t('web.synergy.weak_hint') }}</span>
    </div>
    <p v-if="result.inputs.power_limited" class="notice err">
      {{ t('web.synergy.res_power_warn') }}
    </p>

    <h2 class="sec">{{ t('web.synergy.res_sep') }} · {{ result.pair.drug_a }} × {{ result.pair.drug_b }}</h2>
    <dl class="rows">
      <div class="row"><dt>{{ t('web.synergy.res_s_ab') }}</dt><dd>{{ f(result.separation.s_ab) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_d_aa') }}</dt><dd>{{ f(result.separation.d_aa) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_d_bb') }}</dt><dd>{{ f(result.separation.d_bb) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_d_ab') }}</dt><dd>{{ f(result.separation.d_ab) }}</dd></div>
    </dl>

    <h2 class="sec">{{ t('web.synergy.res_null', { n: result.separation.null.n_perm }) }}</h2>
    <dl class="rows">
      <div class="row"><dt>{{ t('web.synergy.res_null_z') }}</dt><dd>{{ f(result.separation.null.z) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_null_p') }}</dt><dd>{{ f(result.separation.null.p_lower, 3) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_null_mean') }}</dt><dd>{{ f(result.separation.null.mean) }}</dd></div>
      <div class="row"><dt>{{ t('web.synergy.res_null_sd') }}</dt><dd>{{ f(result.separation.null.sd) }}</dd></div>
    </dl>

    <template v-if="result.disease">
      <h2 class="sec">{{ t('web.synergy.res_prox') }}</h2>
      <dl class="rows">
        <div v-for="item in proxRows" :key="item.key" class="row">
          <dt>{{ item.label }}</dt>
          <dd v-if="item.v && item.v.ok">
            z = {{ f(item.v.z, 3) }} · p = {{ f(item.v.p, 3) }} · d_c = {{ f(item.v.d_c, 3) }}
          </dd>
          <dd v-else-if="item.v" class="dim">
            {{ t('web.synergy.res_prox_err', { err: item.v.error ?? '?' }) }}
          </dd>
          <dd v-else class="dim">{{ t('web.synergy.value_na') }}</dd>
        </div>
      </dl>
    </template>

    <p class="meta-line">
      {{ t('web.synergy.res_cov', {
        a: result.inputs.coverage.a_in_lcc, at: result.inputs.coverage.a_total,
        b: result.inputs.coverage.b_in_lcc, bt: result.inputs.coverage.b_total,
        d: result.inputs.coverage.d_in_lcc, dt: result.inputs.coverage.d_total,
        ov: result.inputs.n_overlap_ab }) }}
    </p>
    <p class="meta-line">
      {{ t('web.synergy.res_run', { top: topText(result.meta.run.top_n),
        sort: result.meta.run.sort ?? '—', perm: permText(result.meta.run.proximity_n_perm) }) }}
    </p>
    <p class="meta-line dim">{{ result.meta.computed_at }} · {{ result.meta.engine }}</p>
  </div>
</template>

<style scoped>
.names-form {
  display: flex;
  align-items: flex-end;
  gap: 18px;
}
.names-form .field {
  flex: 1;
}
.disease-block {
  margin-top: 18px;
}
.disease-block .ctl-label {
  display: block;
  margin-bottom: 10px;
}
.param-row {
  display: flex;
  align-items: baseline;
  gap: 12px;
  margin-top: 16px;
  flex-wrap: wrap;
}
.hint {
  font-size: 10.5px;
  letter-spacing: 0.06em;
  color: var(--c-faint);
}
.xlinks {
  margin-top: 20px;
}
.res-kicker {
  margin: 0;
  font-size: 10px;
  letter-spacing: 0.4em;
  color: var(--c-faint);
}
.geometry {
  margin-top: 14px;
  display: flex;
  align-items: baseline;
  gap: 14px;
  flex-wrap: wrap;
}
.geo-label {
  font-size: 17px;
  font-weight: 500;
  letter-spacing: 0.08em;
  color: var(--c-ink);
}
.geo-basis {
  font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace;
  font-size: 11.5px;
  color: var(--c-mid);
}
.weak {
  font-size: 11px;
  letter-spacing: 0.06em;
  color: var(--c-danger);
}
.sec {
  margin: 26px 0 10px;
  font-size: 11px;
  font-weight: 400;
  letter-spacing: 0.3em;
  color: var(--c-sub);
}
.rows {
  margin: 0;
}
.row {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 20px;
  padding: 7px 0;
  border-bottom: 1px solid var(--c-rowline);
}
.row dt {
  font-size: 12px;
  letter-spacing: 0.06em;
  color: var(--c-ink2);
}
.row dd {
  margin: 0;
  font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace;
  font-size: 12px;
  color: var(--c-ink);
}
.meta-line {
  margin: 14px 0 0;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: var(--c-mid);
}
.dim {
  color: var(--c-faint);
}
</style>
