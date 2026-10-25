<script setup lang="ts">
/** 数据管理(阶段 2):背景互作组状态与一键重建(job)+ 分析产物目录(CollectionView)。 */
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type InteractomeStatus } from '../api'
import { useJob } from '../composables/useJob'
import CollectionView from './CollectionView.vue'

const { t } = useI18n()
const itx = ref<InteractomeStatus | null>(null)
const loaded = ref(false)
const notice = ref('')
let noticeTimer = 0

function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 3200)
}

const job = useJob((result) => {
  flash(t('web.data.built_flash', { nodes: result.nodes, edges: result.edges }) +
    (result.cleaned ? t('web.data.built_cleaned', { n: result.cleaned }) : ''))
  void load()
})

async function load(): Promise<void> {
  try {
    itx.value = await api.interactome()
    loaded.value = true
  } catch {
    flash(t('web.data.status_fail'))
  }
}

function build(): void {
  void job.run(api.buildInteractome(), { doneText: t('web.data.build_done') })
}

let pollTimer = 0
onMounted(() => {
  void load()
})
onBeforeUnmount(() => {
  window.clearTimeout(noticeTimer)
  if (!job.running.value) window.clearInterval(pollTimer)
})
</script>

<template>
  <p class="kicker s1">{{ t('web.data.kicker') }}</p>
  <h1 class="title s2">{{ t('web.data.title') }}</h1>
  <p class="sub s3">{{ t('web.data.sub') }}</p>

  <div class="panel s4">
    <div class="controls">
      <span class="ctl-label">{{ t('web.data.itx_label') }}</span>
      <span v-if="loaded && itx?.built" class="count">
        {{ t('web.data.itx_built', { nodes: itx.nodes, edges: itx.edges, version: itx.version,
                                     conf: itx.confidence, at: itx.built_at, age: itx.age ?? '-' }) }}
      </span>
      <span v-else-if="loaded" class="count" style="color:#b3261e">
        {{ t('web.data.itx_missing') }}
      </span>
      <button class="go" :title="t('web.data.rebuild_title')" :disabled="job.running.value" @click="build">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <p v-if="loaded && itx?.built && itx?.stale" class="notice">
      {{ t('web.data.itx_stale', { age: itx.age }) }}
    </p>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ job.text.value || t('web.data.building') }}</p>
    </div>
    <p v-if="job.error.value" class="notice">{{ t('web.data.build_fail', { msg: job.error.value }) }}</p>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>

  <CollectionView kind="downloads" :kicker="t('web.collection.dl_kicker')"
                  :title="t('web.collection.dl_title')" :sub="t('web.collection.dl_sub')" />
</template>
