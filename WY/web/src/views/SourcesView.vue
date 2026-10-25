<script setup lang="ts">
/** 药物/疾病靶点来源管理(阶段 2):清单 + 来源筛选 + 单文件/批量更新(job)+ 打包/删除。
 * 数据源规格唯一定义在 core.sources;更新编排经 /api/sources/{kind}/refresh 以 job 轮询。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type CollectionItem, type CollectionKind } from '../api'
import { useConfirmArm } from '../composables/useConfirmArm'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import { usePager } from '../composables/usePager'
import { useSelection } from '../composables/useSelection'
import { downloadBlob } from '../utils/download'
import { fmtSize } from '../utils/format'

const { t } = useI18n()
const props = defineProps<{ kind: 'drug' | 'disease' }>()

const collectionKind = computed<CollectionKind>(
  () => (props.kind === 'drug' ? 'drug_sources' : 'disease_sources'))

const HEAD = computed(() => ({
  drug: {
    kicker: t('web.sources.drug_kicker'), title: t('web.sources.drug_title'),
    sub: t('web.sources.drug_sub'), zip: 'Drug_sources.zip'
  },
  disease: {
    kicker: t('web.sources.disease_kicker'), title: t('web.sources.disease_title'),
    sub: t('web.sources.disease_sub'), zip: 'Disease_sources.zip'
  }
} as const))
const head = computed(() => HEAD.value[props.kind])

const items = ref<CollectionItem[]>([])
const loaded = ref(false)
const { notice, flash } = useNotice(3200)
const query = ref('')
const { page, paged, nPages, goto } = usePager<CollectionItem>(() => filtered.value)
const { selected, toggle: toggleSel, allSelected, toggleAll, retain } = useSelection(
  () => items.value.map((it) => it.name))

/* 来源筛选 + merged 整合文件开关 */
const srcOn = ref<Record<string, boolean>>({})
const showMerged = ref(true)

const sourceNames = computed(() => {
  const names = new Set<string>()
  for (const it of items.value) {
    if (!it.is_merged && it.source) names.add(it.source)
  }
  return [...names]
})

const filtered = computed<CollectionItem[]>(() => {
  let list = items.value
  if (srcOn.value && Object.keys(srcOn.value).length) {
    list = list.filter((it) => it.is_merged
      ? showMerged.value
      : srcOn.value[it.source ?? ''] !== false)
  }
  const q = query.value.trim().toLowerCase()
  if (q) list = list.filter((it) => it.name.toLowerCase().includes(q))
  return list
})

/* ---------------- 更新(长任务,job 轮询) ---------------- */

const job = useJob((result) => {
  const rows = (result.rows ?? []) as { File: string; Status: string; Message: string }[]
  const failed = rows.filter((r) => r.Status === 'FAILED')
  flash(failed.length
    ? t('web.sources.update_done_fail', { n: failed.length, names: failed.map((r) => r.File).join(', ') })
    : t('web.sources.update_done_ok', { n: rows.length }))
  void load()
})

async function refreshFiles(kind: 'drug' | 'disease', filenames: string[] | null): Promise<void> {
  await job.run(api.refreshSources(kind, filenames), { doneText: t('web.sources.one_ok') })
}

const lastUpdate = ref<{ file: string; status: string; message: string } | null>(null)

async function updateOne(name: string): Promise<void> {
  lastUpdate.value = null
  await refreshFiles(props.kind, [name])
  const rows = (job.result.value?.rows ?? []) as { File: string; Status: string; Message: string }[]
  if (rows.length) lastUpdate.value = { file: rows[0].File, status: rows[0].Status, message: rows[0].Message }
}

/* ---------------- 打包 / 删除 ---------------- */

async function zipSelected(): Promise<void> {
  const names = [...selected.value]
  if (!names.length) return
  try {
    downloadBlob(await api.zip(collectionKind.value, names), head.value.zip)
  } catch {
    flash(t('web.common.zip_fail'))
  }
}

const { armed: confirmArm, tap: confirmTap } = useConfirmArm()

function deleteSelected(): void {
  confirmTap(() => {
    const names = [...selected.value]
    if (!names.length) return
    void (async () => {
      try {
        const r = await api.remove(collectionKind.value, names)
        selected.value = new Set()
        await load()
        flash(r.failed.length
          ? t('web.common.deleted_partial', { n: r.deleted, names: r.failed.join(', ') })
          : t('web.common.deleted', { n: r.deleted }))
      } catch {
        flash(t('web.common.delete_fail'))
      }
    })()
  })
}

/* ---------------- 生命周期 ---------------- */

async function load(): Promise<void> {
  const listing = await api.collection(collectionKind.value)
  items.value = listing.items
  loaded.value = true
  retain(listing.items.map((it) => it.name))
  // 来源筛选项按实际清单初始化(默认全开)
  for (const s of sourceNames.value) {
    if (!(s in srcOn.value)) srcOn.value[s] = true
  }
}

let pollTimer = 0
onMounted(() => {
  void load()
  pollTimer = window.setInterval(() => { if (!job.running.value && !confirmArm.value) void load() }, 15000)
})
onBeforeUnmount(() => {
  window.clearInterval(pollTimer)
})
</script>

<template>
  <p class="kicker s1">{{ head.kicker }}</p>
  <h1 class="title s2">{{ head.title }}</h1>
  <p class="sub s3">{{ head.sub }}</p>

  <div class="panel s4">
    <div class="controls">
      <label class="search">
        <input v-model="query" spellcheck="false" :placeholder="t('web.sources.search_ph')"
               @input="goto(1)" />
        <i class="line"></i>
      </label>
      <span class="count">{{ t('web.common.items_selected', { n: filtered.length, m: selected.size }) }}</span>
      <button class="go" :title="t('web.common.refresh')" :disabled="job.running.value" @click="load">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <div class="filters">
      <span class="ctl-label">{{ t('web.sources.source_label') }}</span>
      <label v-for="s in sourceNames" :key="s" class="chk">
        <input v-model="srcOn[s]" type="checkbox" /> {{ s }}
      </label>
      <label class="chk"><input v-model="showMerged" type="checkbox" /> {{ t('web.sources.merged_label') }}</label>
    </div>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <div v-if="lastUpdate" class="lastupd">
      {{ lastUpdate.file }} — {{ lastUpdate.status === 'OK' ? t('web.sources.one_ok')
        : lastUpdate.status === 'EMPTY' ? t('web.sources.one_empty')
        : t('web.sources.one_fail', { msg: lastUpdate.message }) }}
    </div>

    <div class="bar" :class="{ dim: job.running.value }">
      <button class="bar-btn" :disabled="job.running.value" @click="refreshFiles(props.kind, null)">
        {{ t('web.sources.update_all') }}
      </button>
      <button class="bar-btn" :disabled="job.running.value || !selected.size"
              @click="refreshFiles(props.kind, [...selected].filter((n) => !n.startsWith('merged_')))">
        {{ t('web.sources.update_selected') }}
      </button>
      <button class="bar-btn" :disabled="job.running.value || !selected.size" @click="zipSelected">
        {{ t('web.common.download_selected') }}
      </button>
      <button class="bar-btn danger" :class="{ arm: confirmArm }"
              :disabled="job.running.value || !selected.size" @click="deleteSelected">
        {{ confirmArm ? t('web.common.confirm_delete', { n: selected.size }) : t('web.common.delete_selected') }}
      </button>
    </div>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill" :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ job.text.value || t('web.sources.updating') }}</p>
    </div>
    <p v-if="job.error.value" class="notice">{{ job.error.value }}</p>

    <ul class="files">
      <li v-if="loaded && !filtered.length" class="empty">
        {{ query ? t('web.common.no_match') : t('web.sources.empty') }}
      </li>
      <li v-for="it in paged" :key="it.name" class="file">
        <div class="row">
          <label class="sel">
            <input type="checkbox" :checked="selected.has(it.name)" @change="toggleSel(it.name)" />
          </label>
          <span class="fname" :title="it.name">
            <i v-if="it.is_merged" class="merged-tag">{{ t('web.sources.merged_tag') }}</i>{{ it.name }}
          </span>
          <span class="fmeta">{{ fmtSize(it.size) }} · {{ it.modified }}</span>
          <button class="mini" :title="t('web.sources.refresh_one')" :disabled="job.running.value"
                  @click="updateOne(it.name)">↻</button>
        </div>
      </li>
    </ul>

    <div v-if="nPages > 1" class="pager">
      <button class="pager-btn" :disabled="page <= 1" @click="goto(page - 1)">‹</button>
      <span class="pager-num">{{ page }} / {{ nPages }}</span>
      <button class="pager-btn" :disabled="page >= nPages" @click="goto(page + 1)">›</button>
    </div>
  </div>
</template>

<style scoped>
.lastupd {
  margin: 12px 0 0;
  font-size: 11.5px;
  letter-spacing: 0.06em;
  color: #8a8a88;
}
.merged-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.2em;
  color: #c2c2c2;
  border: 1px solid #e6e6e4;
  padding: 1px 5px;
  margin-right: 8px;
  vertical-align: 1.5px;
}
</style>
