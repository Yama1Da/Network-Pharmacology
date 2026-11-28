<script setup lang="ts">
/** 分析页模块检测折叠区:重跑(Louvain 模块 + 模块富集长表)+ 摘要表与
 * 无信号/失败提示。重跑经 jobs 异步,完成 emit('done', result) 由父级合并回
 * PpiResult.modules(拆分前 onDone 内联合并逻辑);api.rerunModules 自模板
 * 内联移入 script。doneText 走 useJob,完成 flash 在父级。 */
import { ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type PpiResult, type Row } from '../api'
import { useJob } from '../composables/useJob'
import { fmt } from '../utils/format'

const props = defineProps<{
  modules: PpiResult['modules']
  paramsStale: boolean
  directoryRel: string
  prefix: string
}>()
const emit = defineEmits<{
  done: [modules: NonNullable<PpiResult['modules']>]
  /** 重跑忙闲变化:壳层的「开始分析」按钮防重入(重跑中禁点) */
  running: [v: boolean]
}>()
const { t } = useI18n()

const modsJob = useJob((result) => {
  emit('done', toModules(result))
})
watch(() => modsJob.running.value, (v) => emit('running', v))

function rerun(): void {
  if (!props.prefix || props.paramsStale) return
  void modsJob.run(
    api.rerunModules(props.directoryRel, props.prefix),
    { doneText: t('web.analysis.modules_done_text') }
  )
}

/** 重跑结果 → PpiResult.modules 形状(emit 给父级直接合并;字段与后端对齐) */
function toModules(result: Record<string, unknown>): NonNullable<PpiResult['modules']> {
  return {
    summary: (result.summary ?? []) as Row[],
    failures: (result.failures ?? []) as unknown[],
    no_sig: (result.no_sig ?? []) as number[],
    hidden_single: Number(result.hidden_single ?? 0),
    enr_rows: (result.enr_rows ?? []) as Row[],
    enr_total: Number(result.enr_total ?? 0),
    enr_modules: Number(result.enr_modules ?? 0)
  }
}
</script>

<template>
  <section class="sec">
    <details class="fold">
      <summary class="sec-title fold-title">{{ t('web.analysis.sec_modules') }}</summary>
      <button class="bar-btn" :disabled="paramsStale || modsJob.running.value" @click="rerun">
        {{ modsJob.running.value ? t('web.analysis.rerunning') : t('web.analysis.modules_rerun') }}
      </button>
      <p v-if="modsJob.running.value" class="hint-line">{{ modsJob.text.value || t('web.analysis.rerunning') }}</p>
      <p v-if="modsJob.error.value" class="notice err">{{ t('web.analysis.modules_done_text') }}:{{ modsJob.error.value }}</p>
      <p v-if="modules && !modules.summary.length" class="hint-line">
        {{ t('web.analysis.modules_empty') }}
      </p>
      <p v-if="modules?.summary.length" class="hint-line">
        {{ t('web.analysis.modules_summary_title') }}<template v-if="modules.hidden_single">
          · {{ t('web.analysis.modules_hidden', { n: modules.hidden_single }) }}</template>
      </p>
      <div v-if="modules?.summary.length" class="tbl-wrap">
        <table class="tbl">
          <thead><tr><th>module_id</th><th>size</th><th>avg_degree</th><th>genes_preview</th></tr></thead>
          <tbody>
            <tr v-for="r in modules.summary" :key="String(r.module_id)">
              <td>{{ r.module_id }}</td><td>{{ r.size }}</td><td>{{ r.avg_degree }}</td><td>{{ r.genes_preview }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <template v-if="modules?.enr_rows?.length">
        <h4 class="sub-sec">{{ t('web.analysis.mod_enr_title') }}</h4>
        <p class="hint-line">
          {{ t('web.analysis.mod_enr_stats', { total: modules.enr_total }) }}
          <!-- 完整表走底部"下载产物"清单(模块富集完整结果表)与整包 ZIP;此处不再放就地下载(2026-11-28 单一入口定稿) -->
        </p>
        <div class="tbl-wrap">
          <table class="tbl">
            <thead>
              <tr><th>module_id</th><th>source</th><th>term_id</th><th>term_name</th>
                  <th>adj_p</th><th>intersection</th><th>query_size</th></tr>
            </thead>
            <tbody>
              <tr v-for="(r, i) in modules.enr_rows" :key="i">
                <td>{{ r.module_id }}</td><td>{{ r.source }}</td><td>{{ r.term_id }}</td>
                <td>{{ fmt(r.adjusted_p_value, 4) }}</td>
                <td>{{ r.intersection_size }}</td><td>{{ r.query_size }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </template>
      <p v-if="modules?.no_sig?.length" class="hint-line">
        {{ t('web.analysis.modules_no_sig', { names: modules.no_sig.join(', ') }) }}
      </p>
      <p v-if="modules?.failures?.length" class="stale">
        {{ t('web.analysis.modules_failures', { n: modules.failures.length }) }}
      </p>
    </details>
  </section>
</template>
