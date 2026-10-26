<script setup lang="ts">
/** 通用文件集合视图:PPI / 富集 / 收藏 / 数据 四页共用。
 * mode "files" 条目为文件,"dirs" 条目为目录(行内也有独立打包下载);
 * kind=enrichment 附加来源筛选与友好展示名。骨架(顶栏/幕布/页脚)由 App 提供。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type CollectionItem, type CollectionKind } from '../api'
import { useConfirmArm } from '../composables/useConfirmArm'
import { useNotice } from '../composables/useNotice'
import { usePager } from '../composables/usePager'
import { useSelection } from '../composables/useSelection'
import { downloadBlob } from '../utils/download'
import { fmtSize } from '../utils/format'

const { t } = useI18n()
const props = defineProps<{
  kind: CollectionKind
  kicker: string
  title: string
  sub: string
}>()

const items = ref<CollectionItem[]>([])
const mode = ref<'files' | 'dirs'>('files')
const loaded = ref(false)
const { notice, flash } = useNotice(2600)

/* ---------------- 搜索 / 筛选 / 选择 / 分页(客户端) ---------------- */

const query = ref('')

/* 富集页筛选:GSEA 按侧别 + 两来源;无法识别命名的文件不受筛选影响 */
const flt = ref({ gseaDrug: true, gseaDisease: true, gprofiler: true, enrichr: true })

const filtered = computed<CollectionItem[]>(() => {
  let list = items.value
  if (props.kind === 'enrichment') {
    const f = flt.value
    list = list.filter((it) => {
      if (!it.source && !it.side) return true
      if (it.gsea) {
        if (it.side === 'drug') return f.gseaDrug
        if (it.side === 'disease') return f.gseaDisease
        return f.gseaDrug || f.gseaDisease
      }
      if (it.source === 'gProfiler') return f.gprofiler
      if (it.source === 'Enrichr') return f.enrichr
      return true
    })
  }
  const q = query.value.trim().toLowerCase()
  if (q) {
    list = list.filter((it) =>
      it.name.toLowerCase().includes(q) || (it.display ?? '').toLowerCase().includes(q))
  }
  return list
})

const { page, paged, nPages, goto } = usePager<CollectionItem>(() => filtered.value)
const { selected, toggle: toggleSel, allSelected, toggleAll, retain } = useSelection(
  () => items.value.map((it) => it.name))

/* ---------------- 动作 ---------------- */

const { armed: confirmArm, tap: confirmTap } = useConfirmArm()

async function load(): Promise<void> {
  const listing = await api.collection(props.kind)
  items.value = listing.items
  mode.value = listing.mode
  loaded.value = true
  retain(listing.items.map((it) => it.name))
}

async function refresh(): Promise<void> {
  try {
    await load()
    flash(t('web.common.refreshed'))
  } catch {
    flash(t('web.common.list_fail'))
  }
}

async function download(names: string[], single?: string): Promise<void> {
  if (!names.length) return
  try {
    downloadBlob(await api.zip(props.kind, names), single ?? zipName.value)
  } catch {
    flash(t('web.common.zip_fail'))
  }
}

function deleteSelected(): void {
  confirmTap(() => {  // 两段式确认:再次点击才删,4 秒内未点自动撤防
    const names = [...selected.value]
    if (!names.length) return
    void (async () => {
      try {
        const r = await api.remove(props.kind, names)
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

const zipName = computed(() => {
  const map: Record<string, string> = {
    ppi: 'PPI_files.zip', enrichment: 'Enrichment_files.zip',
    stars: 'Starred_results.zip', downloads: 'Data_files.zip'
  }
  return map[props.kind] ?? `${props.kind}.zip`
})

/* ---------------- 生命周期 ---------------- */

let pollTimer = 0
onMounted(() => {
  void load()
  // 目录可能被分析进程写入:静默轻刷新(确认/删除进行中时跳过)
  pollTimer = window.setInterval(() => {
    if (!confirmArm.value) void load()
  }, 15000)
})
onBeforeUnmount(() => {
  window.clearInterval(pollTimer)
})
</script>

<template>
  <p class="kicker s1">{{ kicker }}</p>
  <h1 class="title s2">{{ title }}</h1>
  <p class="sub s3">{{ sub }}</p>

  <div class="panel s4">
    <div class="controls">
      <label class="search">
        <input v-model="query" spellcheck="false" :placeholder="t('web.common.search_ph')"
               @input="goto(1)" />
        <i class="line"></i>
      </label>
      <span class="count">{{ t('web.common.items_selected', { n: filtered.length, m: selected.size }) }}</span>
      <button class="go" :title="t('web.common.refresh')" @click="refresh">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <div v-if="kind === 'enrichment'" class="filters">
      <span class="ctl-label">{{ t('web.collection.flt_label') }}</span>
      <label class="chk"><input v-model="flt.gseaDrug" type="checkbox" /> {{ t('web.collection.flt_gsea_drug') }}</label>
      <label class="chk"><input v-model="flt.gseaDisease" type="checkbox" /> {{ t('web.collection.flt_gsea_disease') }}</label>
      <label class="chk"><input v-model="flt.gprofiler" type="checkbox" /> {{ t('web.collection.flt_gprofiler') }}</label>
      <label class="chk"><input v-model="flt.enrichr" type="checkbox" /> {{ t('web.collection.flt_enrichr') }}</label>
      <span class="ctl-meta">{{ t('web.collection.flt_note') }}</span>
    </div>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <div class="bar" :class="{ dim: !selected.size }">
      <button class="bar-btn" @click="toggleAll">
        {{ allSelected ? t('web.common.unselect_all') : t('web.common.all') }}
      </button>
      <span class="bar-label">{{ selected.size ? t('web.common.selected', { n: selected.size }) : t('web.common.batch_hint') }}</span>
      <button class="bar-btn" :disabled="!selected.size" @click="download([...selected])">{{ t('web.common.download_selected') }}</button>
      <button class="bar-btn danger" :class="{ arm: confirmArm }" :disabled="!selected.size"
              @click="deleteSelected">
        {{ confirmArm ? t('web.common.confirm_delete', { n: selected.size }) : t('web.common.delete_selected') }}
      </button>
    </div>

    <ul class="files">
      <li v-if="loaded && !filtered.length" class="empty">
        {{ query ? t('web.common.no_match') : t('web.common.empty_dir') }}
      </li>
      <li v-for="it in paged" :key="it.name" class="file">
        <div class="row">
          <label class="sel">
            <input type="checkbox" :checked="selected.has(it.name)" @change="toggleSel(it.name)" />
          </label>
          <span class="fname" :title="it.name">
            {{ it.display ?? it.name }}<i v-if="it.display" class="raw">{{ it.name }}</i>
          </span>
          <span class="fmeta">{{ it.size != null ? `${fmtSize(it.size)} · ` : '' }}{{ it.modified }}</span>
          <button class="mini" :title="t('web.common.download')" @click="download([it.name], it.name)">⬇</button>
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
