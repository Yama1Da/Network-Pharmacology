<script setup lang="ts">
/** 分析页 GSEA 折叠区:双侧切换(drug/disease)+ 重跑(参数未过期才可点)+
 * 结果表(NES/adj_p)与产物图。重跑经 jobs 异步,完成 emit('done', result, side)
 * 由父级合并回 PpiResult.gsea(拆分前 onDone 内联合并逻辑);api.rerunGsea
 * 自模板内联移入 script。doneText 走 useJob,完成 flash 在父级。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl, type PpiResult } from '../api'
import { useJob } from '../composables/useJob'
import { fmt } from '../utils/format'

const props = defineProps<{
  gsea: PpiResult['gsea']
  artifacts: Record<string, string>
  paramsStale: boolean
  drug: string
  disease: string
  prefix: string
}>()
const emit = defineEmits<{
  done: [result: Record<string, unknown>, side: 'drug' | 'disease']
  /** 重跑忙闲变化:壳层的「开始分析」按钮防重入(重跑中禁点) */
  running: [v: boolean]
}>()
const { t } = useI18n()

const gseaSide = ref<'drug' | 'disease'>('drug')
const gseaRes = computed(() => props.gsea?.[gseaSide.value] ?? null)
const gseaDone = computed(() => {
  const g = gseaRes.value
  return !!g && !props.paramsStale &&
    !['request_failed', 'circuit_open', 'normalize_empty'].includes(g.reason)
})
const png = computed(() => props.artifacts[`gsea_${gseaSide.value}_png`] ?? null)
watch(() => gseaJob.running.value, (v) => emit('running', v))

const gseaJob = useJob((result) => {
  emit('done', result, gseaSide.value)
})

function rerun(): void {
  if (!props.prefix || gseaDone.value || props.paramsStale) return
  void gseaJob.run(
    api.rerunGsea(gseaSide.value, props.drug, props.disease, props.prefix),
    { doneText: t('web.analysis.gsea_done_text') }
  )
}
</script>

<template>
  <section class="sec">
    <details class="fold">
      <summary class="sec-title fold-title">{{ t('web.analysis.sec_gsea') }}</summary>
      <div class="lv-group">
        <button class="lv" :class="{ on: gseaSide === 'drug' }" @click="gseaSide = 'drug'">{{ t('web.analysis.gsea_side_drug') }}</button>
        <i class="lv-sep">/</i>
        <button class="lv" :class="{ on: gseaSide === 'disease' }" @click="gseaSide = 'disease'">{{ t('web.analysis.gsea_side_disease') }}</button>
      </div>
      <button class="bar-btn" :disabled="gseaDone || paramsStale || gseaJob.running.value" @click="rerun">
        {{ gseaJob.running.value ? t('web.analysis.rerunning') : (gseaDone ? t('web.analysis.gsea_has_result') : t('web.analysis.gsea_rerun')) }}
      </button>
      <p v-if="gseaJob.running.value" class="hint-line">{{ gseaJob.text.value || t('web.analysis.rerunning') }}</p>
      <p v-if="gseaJob.error.value" class="notice err">{{ t('web.analysis.gsea_done_text') }}:{{ gseaJob.error.value }}</p>
      <template v-if="gseaRes">
        <p v-if="gseaRes.reason === 'insufficient_input'" class="stale">{{ t('web.analysis.gsea_insufficient') }}</p>
        <p v-else-if="gseaRes.reason === 'request_failed'" class="stale">{{ t('web.analysis.gsea_req_failed') }}</p>
        <p v-else-if="!gseaRes.rows.length" class="hint-line">{{ t('web.analysis.gsea_empty') }}</p>
        <template v-else>
          <img v-if="png" class="fig-img" :src="artifactsUrl(png)" :alt="`GSEA ${gseaSide}`" />
          <div class="tbl-wrap">
            <table class="tbl">
              <thead><tr><th>term_name</th><th>NES</th><th>adj_p</th></tr></thead>
              <tbody>
                <tr v-for="(r, i) in gseaRes.rows" :key="i">
                  <td>{{ r.Term ?? r.term_name ?? r.Name }}</td>
                  <td>{{ fmt(r.NES ?? r.nes, 3) }}</td>
                  <td>{{ fmt(r['FDR q-val'] ?? r.fdr_q_value ?? r.adj_p ?? r['FDR p-val'], 4) }}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </template>
      </template>
    </details>
  </section>
</template>
