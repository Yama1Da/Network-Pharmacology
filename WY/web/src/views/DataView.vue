<script setup lang="ts">
/** 背景数据库(数据底座,二期随 PK-DB 镜像扩容):背景互作组 + PK-DB 镜像(浓度
 * 可达性数据底座)的状态与一键构建(job)。分析产物目录已拆为独立从属页
 * (#/artifacts),入口在数据网络药理学分析页顶栏(2026-11-17 迁移)。
 * 两面板各自持有独立 notice(2026-11-17 修复:原先共用一个 useNotice,
 * PK 完成提示会在互作组面板里同时冒出来)。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type InteractomeStatus, type PkMirrorStatus } from '../api'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'

const { t } = useI18n()
const itx = ref<InteractomeStatus | null>(null)
const pk = ref<PkMirrorStatus | null>(null)
const loaded = ref(false)
const pkLoaded = ref(false)
/* 查名冲突台账展开(PK 卡):名 → {winner_sid, claimants},按争议名字典序 */
const ledgerOpen = ref(false)
const ledgerRows = computed(() =>
  Object.entries(pk.value?.normalization?.conflicts ?? {})
    .map(([key, v]) => ({ key, winner_sid: v.winner_sid, claimants: v.claimants })))
const { notice: itxNotice, flash: flashItx } = useNotice(3200)
const { notice: pkNotice, flash: flashPk } = useNotice(3200)

const job = useJob((result) => {
  flashItx(t('web.data.built_flash', { nodes: result.nodes, edges: result.edges }) +
    (result.cleaned ? t('web.data.built_cleaned', { n: result.cleaned }) : ''))
  void load()
})

const pkJob = useJob((result) => {
  // 纯续传(0 新建)不再喊"完成:0 研究 / 0 源表",改报"已是最新"
  flashPk(result.studies_built
    ? t('web.data.pk_build_done', { studies: result.studies_built, tsv: result.tsv_files })
    : t('web.data.pk_uptodate', { studies: result.studies_skipped_existing ?? 0 }))
  void load()
})

async function load(): Promise<void> {
  try {
    itx.value = await api.interactome()
    loaded.value = true
    // 在途构建透出(2026-11 A 组):重进页面自动接上进度条,不再靠再点一次按钮
    if (itx.value.building && itx.value.job_id && !job.running.value) {
      void job.resume(itx.value.job_id, { doneText: t('web.data.build_done') })
    }
  } catch {
    flashItx(t('web.data.status_fail'))
  }
  try {
    pk.value = await api.pkMirror()
    pkLoaded.value = true
    if (pk.value.building && pk.value.job_id && !pkJob.running.value) {
      void pkJob.resume(pk.value.job_id, { doneText: t('web.data.build_done') })
    }
  } catch {
    flashPk(t('web.data.pk_status_fail'))
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
      <p v-if="itxNotice" class="notice">{{ itxNotice }}</p>
    </Transition>
  </div>

  <div class="panel s5">
    <div class="controls">
      <span class="ctl-label">{{ t('web.data.pk_label') }}</span>
      <span v-if="pkLoaded && pk?.built" class="count">
        {{ t('web.data.pk_built', { studies: pk.studies, tsv: pk.tsv_files, at: pk.built_at,
                                     age: pk.age ?? '-' }) }}
      </span>
      <span v-if="pkLoaded && pk?.normalization" class="count"
            :title="t('web.data.pk_norm_tip')">
        {{ t('web.data.pk_norm', { v: pk.normalization.pk_norm_version,
                                   c: pk.normalization.conflict_keys ?? 0 }) }}
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

    <!-- 查名冲突台账展开(仅 PK 卡;日志页只有 warning 行,明细在这里) -->
    <div v-if="pkLoaded && (pk?.normalization?.conflict_keys ?? 0) > 0" class="controls">
      <button class="bar-btn" :disabled="pkJob.running.value" @click="ledgerOpen = !ledgerOpen">
        {{ ledgerOpen ? t('web.data.pk_ledger_hide')
           : t('web.data.pk_ledger_btn', { n: pk?.normalization?.conflict_keys ?? 0 }) }}
      </button>
    </div>
    <div v-if="ledgerOpen && ledgerRows.length" class="pk-ledger">
      <div class="ledger-row ledger-head">
        <span>{{ t('web.data.pk_ledger_head_key') }}</span>
        <span>{{ t('web.data.pk_ledger_head_winner') }}</span>
        <span>{{ t('web.data.pk_ledger_head_claimants') }}</span>
      </div>
      <div v-for="row in ledgerRows" :key="row.key" class="ledger-row">
        <span class="mono">{{ row.key }}</span>
        <span class="mono">{{ row.winner_sid }}</span>
        <span class="mono dim">{{ row.claimants.join(', ') }}</span>
      </div>
      <p class="count" style="display:block; color:var(--c-ink3); margin-top:6px">
        {{ t('web.data.pk_ledger_file') }}
      </p>
    </div>
    <p v-if="ledgerOpen && !ledgerRows.length" class="count" style="display:block">
      {{ t('web.data.pk_ledger_empty') }}
    </p>

    <div v-if="pkJob.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((pkJob.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ pkJob.text.value || t('web.data.pk_building') }}</p>
    </div>
    <p v-if="pkJob.error.value" class="notice">{{ t('web.data.build_fail', { msg: pkJob.error.value }) }}</p>
    <Transition name="note">
      <p v-if="pkNotice" class="notice">{{ pkNotice }}</p>
    </Transition>
  </div>
</template>

<style scoped>
.pk-ledger {
  margin-top: 4px;
  border-top: 1px solid var(--c-rowline);
  padding-top: 6px;
  max-height: 280px;
  overflow: auto;
}
.ledger-row {
  display: grid;
  grid-template-columns: minmax(130px, 1.5fr) minmax(110px, 1fr) minmax(120px, 1.2fr);
  gap: 10px;
  padding: 5px 0;
  border-bottom: 1px solid var(--c-rowline);
  font-size: 11.5px;
}
.ledger-head {
  color: var(--c-sub);
  letter-spacing: 0.08em;
  font-size: 10.5px;
}
.mono {
  font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace;
  color: var(--c-ink);
  overflow-wrap: anywhere;
}
.dim {
  color: var(--c-mid);
}
</style>
