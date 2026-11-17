<script setup lang="ts">
/** 数据管理:背景互作组 + PK-DB 镜像(浓度可达性数据底座,二期)的状态与一键构建(job)、分析产物目录。 */
import { onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type InteractomeStatus, type PkMirrorStatus } from '../api'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import CollectionView from './CollectionView.vue'

const { t } = useI18n()
const itx = ref<InteractomeStatus | null>(null)
const pk = ref<PkMirrorStatus | null>(null)
const loaded = ref(false)
const pkLoaded = ref(false)
const { notice, flash } = useNotice(3200)

const job = useJob((result) => {
  flash(t('web.data.built_flash', { nodes: result.nodes, edges: result.edges }) +
    (result.cleaned ? t('web.data.built_cleaned', { n: result.cleaned }) : ''))
  void load()
})

const pkJob = useJob((result) => {
  flash(t('web.data.pk_build_done', { studies: result.studies_built, tsv: result.tsv_files }))
  void load()
})

async function load(): Promise<void> {
  try {
    itx.value = await api.interactome()
    loaded.value = true
  } catch {
    flash(t('web.data.status_fail'))
  }
  try {
    pk.value = await api.pkMirror()
    pkLoaded.value = true
  } catch {
    flash(t('web.data.pk_status_fail'))
  }
}

function build(): void {
  void job.run(api.buildInteractome(), { doneText: t('web.data.build_done') })
}

function buildPk(): void {
  void pkJob.run(api.buildPkMirror(), { doneText: t('web.data.build_done') })
}

onMounted(() => {
  void load()
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
      <span v-else-if="loaded" class="count" style="color:var(--c-danger)">
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

  <div class="panel s5">
    <div class="controls">
      <span class="ctl-label">{{ t('web.data.pk_label') }}</span>
      <span v-if="pkLoaded && pk?.built" class="count">
        {{ t('web.data.pk_built', { studies: pk.studies, tsv: pk.tsv_files, at: pk.built_at,
                                     age: pk.age ?? '-' }) }}
      </span>
      <span v-else-if="pkLoaded" class="count" style="color:var(--c-danger)">
        {{ t('web.data.pk_missing') }}
      </span>
      <button class="go" :title="t('web.data.rebuild_title')" :disabled="pkJob.running.value" @click="buildPk">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <p v-if="pkLoaded && pk?.built && pk?.stale" class="notice">
      {{ t('web.data.pk_stale', { age: pk.age }) }}
    </p>
    <!-- 上游文件集缺失(404 悬空引用):永久状态非故障,灰色弱化,不与失败混淆 -->
    <p v-if="pkLoaded && pk?.built && (pk.missing ?? 0) > 0" class="count"
       style="display:block; color:var(--c-ink3)">
      {{ t('web.data.pk_missing_studies', { n: pk.missing }) }}
    </p>
    <div v-if="pkLoaded && pk?.built && (pk.failures ?? 0) > 0" class="controls">
      <span class="count" style="color:var(--c-danger)">
        {{ t('web.data.pk_failures', { n: pk.failures }) }}
      </span>
      <button class="bar-btn strong" :disabled="pkJob.running.value" @click="buildPk">
        {{ t('web.data.pk_retry_btn') }}
      </button>
    </div>

    <div v-if="pkJob.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((pkJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ pkJob.text.value || t('web.data.pk_building') }}</p>
    </div>
    <p v-if="pkJob.error.value" class="notice">{{ t('web.data.build_fail', { msg: pkJob.error.value }) }}</p>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>

  <CollectionView kind="downloads" :kicker="t('web.collection.dl_kicker')"
                  :title="t('web.collection.dl_title')" :sub="t('web.collection.dl_sub')" />
</template>
