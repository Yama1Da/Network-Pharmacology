<script setup lang="ts">
/** 协同分析结果页(二期,独立管理页):synergy-results 归档清单 · 行级摘要 ·
 * 详情渲染(复用 SynergyResultPanel)· 单跑/批量 ZIP · 删除(两段确认)·
 * 「复现此运行」回填协同预测页表单。数据层在 core.synergy(results_root 等)。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SynergyContract, type SynergyRunRow } from '../api'
import SynergyResultPanel from '../components/SynergyResultPanel.vue'
import { useConfirmArm } from '../composables/useConfirmArm'
import { useNotice } from '../composables/useNotice'
import { usePager } from '../composables/usePager'
import { useSelection } from '../composables/useSelection'
import { setSynergyReplay } from '../store'
import { downloadBlob } from '../utils/download'

const { t } = useI18n()
const runs = ref<SynergyRunRow[]>([])
const runsLoaded = ref(false)
const runsQuery = ref('')
const { notice, flash } = useNotice(4200)
const { selected: runsSel, toggle: runsToggle, allSelected: runsAllSelected,
        toggleAll: runsToggleAll, retain: runsRetain } = useSelection(
  () => runs.value.map((r) => r.run))
const runsDeleteBusy = ref(false)
const runsDeleteErr = ref('')
const { armed: runsConfirm, tap: runsConfirmTap } = useConfirmArm()
const zipBusy = ref(false)
const detail = ref<SynergyContract | null>(null)
const detailRun = ref('')
const detailBusy = ref(false)
const detailErr = ref('')

async function loadRuns(): Promise<void> {
  try {
    runs.value = (await api.synergyResults()).runs
    runsLoaded.value = true
    runsRetain(runs.value.map((r) => r.run))
    // 已删的运行若正展开详情,收起
    if (detailRun.value && !runs.value.some((r) => r.run === detailRun.value)) {
      detailRun.value = ''
      detail.value = null
    }
  } catch (e) {
    flash(t('web.synergy_results.load_fail', { msg: e instanceof Error ? e.message : String(e) }))
  }
}

const runsFiltered = computed(() => {
  const q = runsQuery.value.trim().toLowerCase()
  if (!q) return runs.value
  return runs.value.filter((r) =>
    r.run.toLowerCase().includes(q) ||
    `${r.drug_a} ${r.drug_b}`.toLowerCase().includes(q) ||
    r.disease.toLowerCase().includes(q))
})

const { page: runsPage, paged: runsPaged, nPages: runsNPages, goto: runsGoto } =
  usePager<SynergyRunRow>(() => runsFiltered.value, 4)   // 详情可展开,每页 4 条(2026-11-18 用户定)

async function showDetail(row: SynergyRunRow): Promise<void> {
  if (detailRun.value === row.run) {
    detailRun.value = ''
    detail.value = null
    return
  }
  detailBusy.value = true
  detailErr.value = ''
  try {
    detail.value = await api.synergyResultsDetail(row.run)
    detailRun.value = row.run
  } catch (e) {
    detailErr.value = t('web.synergy_results.detail_fail', {
      msg: e instanceof Error ? e.message : String(e) })
  } finally {
    detailBusy.value = false
  }
}

async function replay(row: SynergyRunRow): Promise<void> {
  // 复现参数取自归档契约的 meta.run/零分布轮数(当时的页级参数),回填即复现
  // (核心确定性已单测验证:同种子同参数 → 同结果)
  let contract: SynergyContract
  try {
    contract = await api.synergyResultsDetail(row.run)
  } catch (e) {
    flash(t('web.synergy_results.detail_fail', {
      msg: e instanceof Error ? e.message : String(e) }))
    return
  }
  setSynergyReplay({
    drug_a: contract.pair.drug_a,
    drug_b: contract.pair.drug_b,
    disease: contract.disease?.name ?? '',
    top_n: contract.meta.run.top_n,
    n_perm: contract.separation.null.n_perm
  })
  flash(t('web.synergy_results.replay_done'))
  location.hash = '#/synergy'
}

async function zipOne(run: string): Promise<void> {
  try {
    downloadBlob(await api.synergyResultsZip([run]), `${run}.zip`)
  } catch {
    flash(t('web.common.zip_fail'))
  }
}

async function zipSelected(): Promise<void> {
  if (!runsSel.value.size) return
  zipBusy.value = true
  try {
    downloadBlob(await api.synergyResultsZip([...runsSel.value]),
                 `synergy_results_${new Date().toISOString().slice(0, 10)}.zip`)
  } catch {
    flash(t('web.common.zip_fail'))
  } finally {
    zipBusy.value = false
  }
}

async function deleteSelected(): Promise<void> {
  runsConfirmTap(() => {
    if (!runsSel.value.size) return
    void (async () => {
      runsDeleteBusy.value = true
      runsDeleteErr.value = ''
      try {
        const r = await api.synergyResultsDelete([...runsSel.value])
        flash(t('web.synergy_results.delete_done', { n: r.deleted.length }))
        runsSel.value = new Set()
        await loadRuns()
      } catch (e) {
        runsDeleteErr.value = e instanceof Error ? e.message : String(e)
      } finally {
        runsDeleteBusy.value = false
      }
    })()
  })
}

const labelOf = (row: SynergyRunRow): string => {
  const m: Record<string, string> = {
    overlap: t('web.synergy.label_overlap'),
    separated: t('web.synergy.label_separated'),
    tie: t('web.synergy.label_tie')
  }
  return m[row.geometry_label] ?? row.geometry_label
}

/* v1.2 E6 分档短词(四档独立,未知≠不可达;旧档/未挂载显示 —);
 * unknown 按真因拆词(2026-11-25):no_exposure=暴露未知,no_activity=活性未知 */
const tierShort = (tier: string | null, reason?: string | null): string => {
  if (tier === 'unknown' && reason === 'no_activity') return t('web.reach.tier_no_activity')
  const m: Record<string, string> = {
    reachable: t('web.synergy.reach_tier_reachable'),
    partial: t('web.synergy.reach_tier_partial'),
    unknown: t('web.synergy.reach_tier_unknown'),
    unreachable: t('web.synergy.reach_tier_unreachable')
  }
  return (tier && m[tier]) || t('web.synergy.value_na')
}

onMounted(() => void loadRuns())
</script>

<template>
  <p class="kicker s1">{{ t('web.synergy_results.kicker') }}</p>
  <h1 class="title s2">{{ t('web.synergy_results.title') }}</h1>
  <p class="sub s3">{{ t('web.synergy_results.sub') }}</p>

  <div class="panel wide s4">
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>

    <p v-if="runsLoaded && !runs.length" class="empty">{{ t('web.synergy_results.empty') }}</p>

    <div class="controls">
      <label class="search">
        <input v-model="runsQuery" spellcheck="false" :placeholder="t('web.synergy_results.search_ph')"
               @input="runsGoto(1)" />
        <i class="line"></i>
      </label>
      <span class="count">{{ t('web.common.selected', { n: runsSel.size }) }}</span>
      <button class="go" :title="t('web.common.refresh')" @click="loadRuns">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <div class="bar" :class="{ dim: !runsSel.size }">
      <button class="bar-btn" @click="runsToggleAll">{{ runsAllSelected ? t('web.common.unselect_all') : t('web.common.all') }}</button>
      <span class="bar-label">{{ runsSel.size ? t('web.common.selected', { n: runsSel.size }) : t('web.common.batch_hint') }}</span>
      <button class="bar-btn" :disabled="!runsSel.size || zipBusy" @click="zipSelected">
        {{ t('web.common.download_selected') }}
      </button>
      <button class="bar-btn danger" :class="{ arm: runsConfirm }" :disabled="!runsSel.size || runsDeleteBusy"
              @click="deleteSelected">
        {{ runsConfirm ? t('web.common.confirm_delete', { n: runsSel.size }) : t('web.common.delete_selected') }}
      </button>
    </div>
    <p v-if="runsDeleteErr" class="notice err">{{ runsDeleteErr }}</p>

    <ul class="files">
      <li v-if="runsLoaded && !runsFiltered.length" class="empty">{{ t('web.common.no_match') }}</li>
      <li v-for="r in runsPaged" :key="r.run" class="file">
        <div class="row">
          <label class="sel">
            <input type="checkbox" :checked="runsSel.has(r.run)" @change="runsToggle(r.run)" />
          </label>
          <span class="fname" :title="r.run">{{ r.drug_a }} × {{ r.drug_b
          }}<template v-if="r.disease"> · {{ r.disease }}</template><i v-if="!r.ok" class="bad-tag">{{ r.reason }}</i></span>
          <span class="fmeta">
            {{ labelOf(r) }} · s_AB {{ r.s_ab ?? '—' }}
            <template v-if="r.r_tier_a || r.r_tier_b">
              · {{ t('web.synergy_results.reach_tag', { a: tierShort(r.r_tier_a, r.r_unknown_a), b: tierShort(r.r_tier_b, r.r_unknown_b) }) }}
            </template>
            <i v-if="r.r_high_conf" class="conf-tag">{{ t('web.synergy_results.high_conf_tag') }}</i>
            <i v-if="r.r_sort_key_note === 'z_null'" class="limit-tag">{{ t('web.synergy.reach_sort_null') }}</i>
            <i v-if="r.power_limited" class="limit-tag">{{ t('web.synergy_results.power_tag') }}</i>
            · {{ r.mtime }}
          </span>
        </div>
        <div class="row-actions">
          <button class="bar-btn" type="button" :disabled="!r.ok" @click="showDetail(r)">
            {{ detailRun === r.run ? t('web.synergy_results.detail_hide') : t('web.synergy_results.detail') }}
          </button>
          <button class="bar-btn" type="button" :disabled="!r.ok" @click="replay(r)">
            {{ t('web.synergy_results.replay') }}
          </button>
          <button class="bar-btn" type="button" @click="zipOne(r.run)">ZIP</button>
        </div>
        <div v-if="detailRun === r.run && detail" class="detail">
          <SynergyResultPanel :result="detail" />
        </div>
        <p v-if="detailRun === r.run && detailErr" class="notice err">{{ detailErr }}</p>
      </li>
    </ul>

    <div v-if="runsNPages > 1" class="pager">
      <button class="pager-btn" :disabled="runsPage <= 1" @click="runsGoto(runsPage - 1)">‹</button>
      <span class="pager-num">{{ runsPage }} / {{ runsNPages }}</span>
      <button class="pager-btn" :disabled="runsPage >= runsNPages" @click="runsGoto(runsPage + 1)">›</button>
    </div>
  </div>
</template>

<style scoped>
.bad-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.15em;
  color: var(--c-danger);
  border: 1px solid var(--c-danger);
  padding: 0 5px;
  margin-left: 10px;
  vertical-align: 1.5px;
}
.limit-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.15em;
  color: var(--c-sub);
  border: 1px solid var(--c-faint);
  padding: 0 5px;
  margin-left: 8px;
  vertical-align: 1.5px;
}
.conf-tag {
  font-style: normal;
  font-size: 10px;
  letter-spacing: 0.15em;
  color: var(--c-accent);
  border: 1px solid currentColor;
  padding: 0 5px;
  margin-left: 8px;
  vertical-align: 1.5px;
}
.row-actions {
  display: flex;
  align-items: baseline;
  gap: 14px;
  padding: 6px 0 2px;
  border-bottom: 1px solid var(--c-rowline);
}
.detail {
  margin-top: 14px;
  padding: 4px 2px 10px;
  border-left: 2px solid var(--c-line);
  padding-left: 16px;
}
</style>
