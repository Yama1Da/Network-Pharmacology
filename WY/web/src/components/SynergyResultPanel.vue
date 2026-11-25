<script setup lang="ts">
/** 协同结果渲染面板(契约 synergy.v1.2):SynergyView 与协同分析结果页两处复用。
 * 措辞红线(§1.3):三枚举只述几何事实,禁译"协同/拮抗";|z_sep|<1 附弱信号提示。
 * v1.2 E6:可达性四档为药理学事实,同样不做协同/拮抗断言;"暴露未知"独立态,
 * 不得与"不可达"合并显示(§2 HomeView 行);fub 保守上界随行标注"可能高估"。
 * 深色适配:色值全走 var(--c-*),无硬编码。 */
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import type { SynergyContract, SynergyProx, SynergyReachDrug } from '../api'
import { fmt } from '../utils/format'

const props = defineProps<{ result: SynergyContract }>()
const { t } = useI18n()

const labelText = computed(() => {
  const m: Record<string, string> = {
    overlap: t('web.synergy.label_overlap'),
    separated: t('web.synergy.label_separated'),
    tie: t('web.synergy.label_tie')
  }
  return m[props.result.geometry_label] ?? props.result.geometry_label
})

/* 弱信号(2026-11-18 实证:宽靶点集/噪声线附近标签随符号翻转):零分布退化或 |z_sep|<1 */
const weakHint = computed(() => {
  const s = props.result.separation
  return s.null.z === null || Math.abs(s.null.z) < 1
})

/* 同药组合:结果为平凡重叠(恒等式 s(A,A)=−d_AA/n),仅供自检基线,不作结论 */
const samePair = computed(() => {
  const a = props.result.pair.drug_a.trim().toLowerCase()
  const b = props.result.pair.drug_b.trim().toLowerCase()
  return a !== '' && a === b
})

const f = (x: number | null | undefined, digits = 4): string =>
  x === null || x === undefined ? t('web.synergy.value_na') : fmt(x, digits)

const topText = (n: number | null): string =>
  n === null ? t('web.synergy.param_n_full') : String(n)

const permText = (n: number | null): string =>
  n === null ? t('web.synergy.value_na') : String(n)

/* proximity 三行渲染:ok=true 显示数值;否则透出后端错误令牌(不伪装成数值) */
const proxRows = computed(() => {
  const p = props.result.proximity
  if (!p || !props.result.disease) return []
  const rows: { key: string; label: string; v: SynergyProx | null }[] = [
    { key: 'z_a', label: t('web.synergy.res_prox_z_a'), v: p.z_a },
    { key: 'z_b', label: t('web.synergy.res_prox_z_b'), v: p.z_b },
    { key: 'z_ab', label: t('web.synergy.res_prox_z_ab'), v: p.z_ab }
  ]
  return rows
})

/* ---- v1.2 E6 可达性块(null = 未挂载/缓存未就绪/旧档,整节不渲染) ---- */
const reach = computed(() => props.result.reachability ?? null)

/* unknown 按真因拆词(2026-11-25):no_activity=活性未知,其余 unknown=暴露未知 */
const tierText = (tier: string | null | undefined, reason?: string | null): string => {
  if (tier === 'unknown' && reason === 'no_activity') return t('web.reach.tier_no_activity')
  const m: Record<string, string> = {
    reachable: t('web.synergy.reach_tier_reachable'),
    partial: t('web.synergy.reach_tier_partial'),
    unknown: t('web.synergy.reach_tier_unknown'),
    unreachable: t('web.synergy.reach_tier_unreachable')
  }
  return (tier && m[tier]) || tier || t('web.synergy.value_na')
}

/* 单药摘要行:分档 · 可算比例 · 游离 Cmax 中位(µM);空值如实显示 — */
const drugLine = (d: SynergyReachDrug | undefined): string => {
  if (!d) return t('web.synergy.value_na')
  const parts = [tierText(d.tier, d.unknown_reason),
                 t('web.synergy.reach_computed', { n: d.n_computed, total: d.n_targets })]
  const cmax = d.exposure?.cmax_molar_M
  if (cmax !== null && cmax !== undefined) {
    parts.push(`Cmax_free ${fmt(cmax * 1e6, 2)} µM`)
  }
  return parts.join(' · ')
}

/* 兜底标记(契约 v1.1):exposure.synonym_fallback_sid 非空 = 数据取自同义兄弟
 * 节点(主 sid 无暴露记录时的台账兜底)——JSON 已存,界面必须可见才算"显式" */
const fallbackSids = computed(() => ({
  a: reach.value?.drug_a?.exposure?.synonym_fallback_sid ?? null,
  b: reach.value?.drug_b?.exposure?.synonym_fallback_sid ?? null
}))
</script>

<template>
  <p class="res-kicker">{{ t('web.synergy.res_title') }} · {{ result.schema }}</p>

  <div class="geometry">
    <span class="geo-label">{{ labelText }}</span>
    <span class="geo-basis">{{ result.label_basis }}</span>
    <span v-if="weakHint" class="weak">{{ t('web.synergy.weak_hint') }}</span>
  </div>
  <p v-if="samePair" class="notice err">{{ t('web.synergy.same_pair_hint') }}</p>
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

  <template v-if="reach">
    <h2 class="sec">{{ t('web.synergy.res_reach') }}</h2>
    <dl class="rows">
      <div class="row">
        <dt>{{ result.pair.drug_a }}</dt>
        <dd>{{ drugLine(reach.drug_a) }}<i v-if="reach.drug_a?.exposure?.fub_imputed"
               class="reach-warn">{{ t('web.synergy.reach_fub_imputed') }}</i><i
               v-if="fallbackSids.a" class="reach-warn">{{ t('web.synergy.reach_synonym_fb', { sid: fallbackSids.a }) }}</i></dd>
      </div>
      <div class="row">
        <dt>{{ result.pair.drug_b }}</dt>
        <dd>{{ drugLine(reach.drug_b) }}<i v-if="reach.drug_b?.exposure?.fub_imputed"
               class="reach-warn">{{ t('web.synergy.reach_fub_imputed') }}</i><i
               v-if="fallbackSids.b" class="reach-warn">{{ t('web.synergy.reach_synonym_fb', { sid: fallbackSids.b }) }}</i></dd>
      </div>
      <div class="row">
        <dt>{{ t('web.synergy.reach_high_conf_label') }}</dt>
        <dd>{{ reach.high_confidence ? t('web.synergy.reach_conf_yes') : t('web.synergy.reach_conf_no') }}</dd>
      </div>
      <div class="row">
        <dt>{{ t('web.synergy.reach_sort_key') }}</dt>
        <dd v-if="reach.sort_key !== null && reach.sort_key !== undefined">{{ fmt(reach.sort_key, 4) }}</dd>
        <dd v-else class="dim">{{ t('web.synergy.reach_sort_null') }}</dd>
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
</template>

<style scoped>
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
.reach-warn {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.1em;
  color: var(--c-mid);
  border: 1px solid var(--c-faint);
  padding: 0 4px;
  margin-left: 8px;
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
