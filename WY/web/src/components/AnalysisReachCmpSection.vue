<script setup lang="ts">
/** 分析页暴露约束前后排序对比(W7 招牌图,2026-11-27):门槛 max(5, ceil(5%×N)),
 * 过门槛看 artifacts 的 reach_compare_png;低于门槛只出结论行(方案 §3.2)。
 * 纯展示,gate 摘要与产物表由父级下发。 */
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { artifactsUrl } from '../api'

const props = defineProps<{
  reachCmp: { gate: boolean; judged: number; total: number }
  artifacts: Record<string, string>
}>()
const { t } = useI18n()

const png = computed(() => props.artifacts['reach_compare_png'] ?? null)
</script>

<template>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.analysis.reach_cmp_title') }}</h3>
    <template v-if="png">
      <p class="hint-line">{{ t('web.analysis.reach_cmp_hint') }}</p>
      <img class="fig-img" :src="artifactsUrl(png)" alt="reachability-weighted rank comparison" />
    </template>
    <p v-else class="hint-line">
      {{ t('web.analysis.reach_cmp_gate', { judged: reachCmp.judged, total: reachCmp.total }) }}
    </p>
  </section>
</template>
