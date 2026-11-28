<script setup lang="ts">
/** 分析页邻近度检验展示:三元指标(Dc/z/p)+ 排除枢纽的排除版指标 + 覆盖与
 * 缓存命中说明。纯只读投影,数据由父级从 PpiResult.proximity 下发;names 模式
 * 的渲染条件由壳层控制。not_built / 低覆盖 / 失败三态独立提示。 */
import { useI18n } from 'vue-i18n'
import type { ProximityResult } from '../api'
import { fmt } from '../utils/format'

defineProps<{ prox: ProximityResult | null }>()
const { t } = useI18n()
</script>

<template>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.analysis.sec_prox') }}</h3>
    <p class="hint-line">
      {{ t('web.analysis.prox_setup', { n: prox?.n_perm, seed: prox?.seed ?? '-' }) }}
    </p>
    <template v-if="prox?.error === 'not_built'">
      <p class="stale">{{ t('web.analysis.prox_not_built') }}</p>
    </template>
    <template v-else-if="prox?.error &&
                         ['low_coverage_drug', 'low_coverage_disease', 'empty_after_filter'].includes(prox.error)">
      <p class="stale">
        {{ t('web.analysis.prox_low_cov', { din: prox.coverage?.drug_in_lcc ?? 0, dtotal: prox.coverage?.drug_total ?? 0,
                                             iin: prox.coverage?.disease_in_lcc ?? 0, itotal: prox.coverage?.disease_total ?? 0 }) }}
      </p>
    </template>
    <template v-else-if="prox?.error">
      <p class="stale">{{ t('web.analysis.prox_fail', { err: prox.error }) }}</p>
    </template>
    <template v-else-if="prox?.ok">
      <div class="metrics three">
        <div class="metric"><b>{{ fmt(prox.d_c) }}</b><span>{{ t('web.analysis.m_dc') }}</span></div>
        <div class="metric"><b>{{ fmt(prox.z) }}</b><span>{{ t('web.analysis.m_z') }}</span></div>
        <div class="metric"><b>{{ fmt(prox.p_empirical, 4) }}</b><span>{{ t('web.analysis.m_p') }}</span></div>
      </div>
      <p class="hint-line" :class="{ good: typeof prox.z === 'number' && prox.z < 0 }">
        {{ typeof prox.z === 'number'
          ? (prox.z < 0 ? t('web.analysis.prox_support') : t('web.analysis.prox_unsupported'))
          : t('web.analysis.prox_degenerate') }}
      </p>
      <p class="hint-line">
        {{ t('web.analysis.prox_lcc', { nodes: prox.lcc_nodes, din: prox.coverage?.drug_in_lcc ?? 0,
                                        dtotal: prox.coverage?.drug_total ?? 0,
                                        iin: prox.coverage?.disease_in_lcc ?? 0,
                                        itotal: prox.coverage?.disease_total ?? 0 }) }}
        {{ prox.cached ? t('web.analysis.cached_hit') : '' }}
      </p>
      <template v-if="prox.d_c_excl != null">
        <p class="hint-line">{{ t('web.analysis.prox_excl', { n: prox.n_excl ?? 0 }) }}</p>
        <div class="metrics three">
          <div class="metric"><b>{{ fmt(prox.d_c_excl) }}</b><span>{{ t('web.analysis.m_dc_excl') }}</span></div>
          <div class="metric"><b>{{ fmt(prox.z_excl) }}</b><span>{{ t('web.analysis.m_z_excl') }}</span></div>
          <div class="metric"><b>{{ fmt(prox.p_excl, 4) }}</b><span>{{ t('web.analysis.m_p_excl') }}</span></div>
        </div>
      </template>
    </template>
  </section>
</template>
