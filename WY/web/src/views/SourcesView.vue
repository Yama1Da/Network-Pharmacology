<script setup lang="ts">
/** 药物/疾病靶点来源管理(阶段 2):清单 + 来源筛选 + 单文件/批量更新(job)+ 打包/删除。
 * 数据源规格唯一定义在 core.sources;更新编排经 /api/sources/{kind}/refresh 以 job 轮询。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { api, type CollectionItem, type CollectionKind } from '../api'
import { useJob } from '../composables/useJob'

const props = defineProps<{ kind: 'drug' | 'disease' }>()

const collectionKind = computed<CollectionKind>(
  () => (props.kind === 'drug' ? 'drug_sources' : 'disease_sources'))

const HEAD = {
  drug: {
    kicker: 'DRUG TARGET SOURCES', title: '药物靶点来源',
    sub: 'ChEMBL + DGIdb · 更新后自动重建整合表', zip: 'Drug_sources.zip'
  },
  disease: {
    kicker: 'DISEASE TARGET SOURCES', title: '疾病靶点来源',
    sub: 'Open Targets + DisGeNET · 更新后自动重建整合表', zip: 'Disease_sources.zip'
  }
} as const
const head = computed(() => HEAD[props.kind])

const items = ref<CollectionItem[]>([])
const loaded = ref(false)
const query = ref('')
const selected = ref<Set<string>>(new Set())
const page = ref(1)
const PAGE_SIZE = 10

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

const paged = computed(() =>
  filtered.value.slice((page.value - 1) * PAGE_SIZE, page.value * PAGE_SIZE))
const nPages = computed(() => Math.max(1, Math.ceil(filtered.value.length / PAGE_SIZE)))
function goto(p: number): void {
  page.value = Math.min(Math.max(1, p), nPages.value)
}

function toggleSel(name: string): void {
  const next = new Set(selected.value)
  if (next.has(name)) next.delete(name)
  else next.add(name)
  selected.value = next
}
const allSelected = computed(
  () => items.value.length > 0 && items.value.every((it) => selected.value.has(it.name)))
function toggleAll(): void {
  selected.value = allSelected.value ? new Set() : new Set(items.value.map((it) => it.name))
}

/* ---------------- 更新(长任务,job 轮询) ---------------- */

const job = useJob((result) => {
  const rows = (result.rows ?? []) as { File: string; Status: string; Message: string }[]
  const failed = rows.filter((r) => r.Status === 'FAILED')
  flash(failed.length
    ? `更新完成,${failed.length} 个文件失败:${failed.map((r) => r.File).join(', ')}`
    : `更新完成:${rows.length} 个文件全部成功`)
  void load()
})

const notice = ref('')
let noticeTimer = 0
function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 3200)
}

async function refreshFiles(kind: 'drug' | 'disease', filenames: string[] | null): Promise<void> {
  await job.run(api.refreshSources(kind, filenames), { doneText: '更新完成' })
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
    const blob = await api.zip(collectionKind.value, names)
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = head.value.zip
    a.click()
    URL.revokeObjectURL(url)
  } catch {
    flash('打包下载失败')
  }
}

const confirmArm = ref(false)
let confirmTimer = 0
function deleteSelected(): void {
  if (!confirmArm.value) {
    confirmArm.value = true
    window.clearTimeout(confirmTimer)
    confirmTimer = window.setTimeout(() => (confirmArm.value = false), 4000)
    return
  }
  confirmArm.value = false
  const names = [...selected.value]
  if (!names.length) return
  void (async () => {
    try {
      const r = await api.remove(collectionKind.value, names)
      selected.value = new Set()
      await load()
      flash(r.failed.length ? `删除 ${r.deleted} 项,失败:${r.failed.join(', ')}` : `已删除 ${r.deleted} 项`)
    } catch {
      flash('删除失败')
    }
  })()
}

/* ---------------- 生命周期 ---------------- */

async function load(): Promise<void> {
  const listing = await api.collection(collectionKind.value)
  items.value = listing.items
  loaded.value = true
  const names = new Set(listing.items.map((it) => it.name))
  selected.value = new Set([...selected.value].filter((n) => names.has(n)))
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
  window.clearTimeout(noticeTimer)
  window.clearTimeout(confirmTimer)
  window.clearInterval(pollTimer)
})

const fmtSize = (n: number): string =>
  n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(1)} KB` : `${(n / 1048576).toFixed(1)} MB`
</script>

<template>
  <p class="kicker s1">{{ head.kicker }}</p>
  <h1 class="title s2">{{ head.title }}</h1>
  <p class="sub s3">{{ head.sub }}</p>

  <div class="panel s4">
    <div class="controls">
      <label class="search">
        <input v-model="query" spellcheck="false" placeholder="搜索文件名…" @input="goto(1)" />
        <i class="line"></i>
      </label>
      <span class="count">{{ filtered.length }} 项 · 已选 {{ selected.size }}</span>
      <button class="go" title="刷新" :disabled="job.running.value" @click="load">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <div class="filters">
      <span class="ctl-label">来源</span>
      <label v-for="s in sourceNames" :key="s" class="chk">
        <input v-model="srcOn[s]" type="checkbox" /> {{ s }}
      </label>
      <label class="chk"><input v-model="showMerged" type="checkbox" /> 整合表 (merged)</label>
    </div>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <div v-if="lastUpdate" class="lastupd">
      {{ lastUpdate.file }} — {{ lastUpdate.status === 'OK' ? '更新成功'
        : lastUpdate.status === 'EMPTY' ? '抓取成功但 0 靶点(检查查询物)'
        : `失败:${lastUpdate.message}` }}
    </div>

    <div class="bar" :class="{ dim: job.running.value }">
      <button class="bar-btn" :disabled="job.running.value" @click="refreshFiles(props.kind, null)">
        全部更新
      </button>
      <button class="bar-btn" :disabled="job.running.value || !selected.size"
              @click="refreshFiles(props.kind, [...selected].filter((n) => !n.startsWith('merged_')))">
        更新选中
      </button>
      <button class="bar-btn" :disabled="job.running.value || !selected.size" @click="zipSelected">
        下载选中 (ZIP)
      </button>
      <button class="bar-btn danger" :class="{ arm: confirmArm }"
              :disabled="job.running.value || !selected.size" @click="deleteSelected">
        {{ confirmArm ? `确认删除 ${selected.size} 项?` : '删除选中' }}
      </button>
    </div>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill" :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ job.text.value || '更新中…' }}</p>
    </div>
    <p v-if="job.error.value" class="notice">{{ job.error.value }}</p>

    <ul class="files">
      <li v-if="loaded && !filtered.length" class="empty">
        {{ query ? '搜索无匹配。' : '尚无来源文件 —— 完成一次分析后自动生成。' }}
      </li>
      <li v-for="it in paged" :key="it.name" class="file">
        <div class="row">
          <label class="sel">
            <input type="checkbox" :checked="selected.has(it.name)" @change="toggleSel(it.name)" />
          </label>
          <span class="fname" :title="it.name">
            <i v-if="it.is_merged" class="merged-tag">整合</i>{{ it.name }}
          </span>
          <span class="fmeta">{{ fmtSize(it.size) }} · {{ it.modified }}</span>
          <button class="mini" title="立即更新" :disabled="job.running.value"
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
