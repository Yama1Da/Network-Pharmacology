<script setup lang="ts">
/** 分析页 RWR 双侧结果(疾病种子 / 药物种子折叠区):迭代元信息、度匹配校正
 * 富集检验、非种子 Top 表 + 种子清单折叠。纯只读投影,数据由父级从
 * PpiResult.rwr 下发;names 模式的渲染条件由壳层控制。 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import type { RwrRow, RwrSideResult } from '../api'
import { fmt } from '../utils/format'

const props = defineProps<{
  rwr: Record<string, RwrSideResult> | null
}>()
const { t } = useI18n()

const rwrSide = ref<'disease_seeds' | 'drug_seeds'>('disease_seeds')
const rwrRes = computed(() => props.rwr?.[rwrSide.value] ?? null)
const rwrTop = computed<(RwrRow & { is_seed: boolean })[]>(() => {
  const res = rwrRes.value
  if (!res?.ok) return []
  const seeds = (res.seeds ?? []).map((r) => ({ ...r, is_seed: true }))
  const non = (res.top_nonseed ?? []).map((r) => ({ ...r, is_seed: false }))
  return [...seeds, ...non].sort((a, b) => b.score - a.score)
})
const rwrNonSeed = computed(() => rwrTop.value.filter((r) => !r.is_seed))
const rwrSeeds = computed(() => rwrTop.value.filter((r) => r.is_seed))
</script>

<template>
  <section class="sec">
    <details class="fold">
      <summary class="sec-title fold-title">{{ t('web.analysis.sec_rwr') }}</summary>
      <p class="hint-line">{{ t('web.analysis.rwr_hint') }}</p>
      <div class="lv-group">
        <button class="lv" :class="{ on: rwrSide === 'disease_seeds' }"
                @click="rwrSide = 'disease_seeds'">{{ t('web.analysis.rwr_side_disease') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: rwrSide === 'drug_seeds' }"
                @click="rwrSide = 'drug_seeds'">{{ t('web.analysis.rwr_side_drug') }}</button>
      </div>
      <template v-if="rwrRes">
        <p v-if="rwrRes.error === 'not_built'" class="stale">{{ t('web.analysis.rwr_not_built') }}</p>
        <p v-else-if="rwrRes.error" class="stale">{{ t('web.analysis.rwr_fail', { err: rwrRes.error }) }}</p>
        <template v-else-if="rwrRes.ok">
          <p class="hint-line">
            {{ t('web.analysis.rwr_meta', { n: rwrRes.n_iter,
                                            conv: rwrRes.converged ? t('web.analysis.converged') : t('web.analysis.not_converged'),
                                            in: rwrRes.coverage?.seed_in_lcc ?? 0, total: rwrRes.coverage?.seed_total ?? 0 }) }}
            {{ rwrRes.cached ? t('web.analysis.cached_hit') : '' }}
          </p>
          <p v-if="rwrRes.enrichment && rwrRes.enrichment.z != null" class="hint-line"
             :class="{ good: rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05 }">
            {{ t('web.analysis.rwr_enrich', {
              n: rwrRes.enrichment.n_annotated, z: fmt(rwrRes.enrichment.z),
              p: fmt(rwrRes.enrichment.p, 4),
              sig: rwrRes.enrichment.p != null && rwrRes.enrichment.p < 0.05
                ? t('web.analysis.sig_yes') : t('web.analysis.sig_no') }) }}
          </p>
          <div v-if="rwrNonSeed.length" class="tbl-wrap">
            <table class="tbl">
              <thead><tr><th>gene</th><th>score</th><th v-if="'annotated' in rwrNonSeed[0]">{{ t('web.analysis.th_annot') }}</th></tr></thead>
              <tbody>
                <tr v-for="r in rwrNonSeed" :key="String(r.gene)">
                  <td>{{ r.gene }}</td><td>{{ fmt(r.score, 6) }}</td>
                  <td v-if="'annotated' in r">{{ r.annotated ? '✓' : '' }}</td>
                </tr>
              </tbody>
            </table>
          </div>
          <p v-else class="hint-line">{{ t('web.analysis.rwr_empty') }}</p>
          <details v-if="rwrSeeds.length" class="fold">
            <summary class="hint-line">{{ t('web.analysis.seeds_fold', { n: rwrSeeds.length }) }}</summary>
            <div class="tbl-wrap">
              <table class="tbl">
                <thead><tr><th>gene</th><th>score</th></tr></thead>
                <tbody><tr v-for="r in rwrSeeds" :key="String(r.gene)"><td>{{ r.gene }}</td><td>{{ fmt(r.score, 6) }}</td></tr></tbody>
              </table>
            </div>
          </details>
        </template>
      </template>
    </details>
  </section>
</template>
