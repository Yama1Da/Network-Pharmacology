import { computed, ref, type Ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type Row } from '../api'
import { freshBadge, pruneBadges, putBadge } from '../utils/probeBadgeCache'

/** 检索表体检徽章批检(MlView ① 区):行级预检状态存 utils/probeBadgeCache
 * (模块级,切页不丢;判定是 GEO 归档不可变属性,缓存安全;本地盘上状态变化
 * 由 overviewJob/forceJob 的 onDone 钩子精准失效,TTL 1h 兜底)+ 勾选集 +
 * 批检忙碌。flash 为壳层 useNotice 的页顶提示。在 setup 上下文中调用。 */
export function useMlProbeBadges(searchResults: Ref<Row[]>, flash: (msg: string) => void) {
  const { t } = useI18n()
  const selectedRows = ref<Set<string>>(new Set())
  const probeBatchBusy = ref(false)
  const selectedCount = computed(() => selectedRows.value.size)
  const allSelected = computed(() =>
    searchResults.value.length > 0 &&
    searchResults.value.every((r) => selectedRows.value.has(rowGse(r))))

  function rowGse(r: Row): string {
    return String(r.gse ?? '').trim().toUpperCase()
  }

  function toggleRow(g: string): void {
    const s = new Set(selectedRows.value)
    if (s.has(g)) s.delete(g)
    else s.add(g)
    selectedRows.value = s
  }

  function toggleAll(): void {
    selectedRows.value = allSelected.value
      ? new Set()
      : new Set(searchResults.value.map(rowGse))
  }

  /* 体检后的 symbol_level 实测结论:i18n 键 + 语义色(空串=维持启发式原文) */
  function symbolVerdict(g: string): { key: string; tone: 'ok' | 'bad' | 'unknown' } | null {
    const p = freshBadge(g)
    if (!p) return null
    if (p.status === 'err') return { key: 'web.ml.sym_unknown', tone: 'unknown' }
    const a = p.assess
    if (!a) return { key: 'web.ml.sym_unknown', tone: 'unknown' }
    if (a.matrix_empty && a.counts_available) return { key: 'web.ml.sym_ok_counts', tone: 'ok' }
    if (a.matrix_empty) return { key: 'web.ml.sym_bad_no_counts', tone: 'bad' }
    if (a.symbol_level) return { key: 'web.ml.sym_ok_gene', tone: 'ok' }
    if (a.annot === 'missing') return { key: 'web.ml.sym_bad_annot', tone: 'bad' }
    return { key: 'web.ml.sym_unknown', tone: 'unknown' }
  }

  function dotClass(g: string): string {
    // 批检进行中:已勾选且未出结果的行亮琥珀脉冲(结果由同步端点一次性返回,
    // 期间全部选中行保持 busy 态,到齐后统一翻转为终态)
    if (!freshBadge(g)) {
      return probeBatchBusy.value && selectedRows.value.has(g) ? 'busy' : 'unknown'
    }
    const v = symbolVerdict(g)!
    return v.tone === 'ok' ? 'ok' : v.tone === 'bad' ? 'bad' : 'unknown'
  }

  /** 检索表行的 platform/type(GSE → 平台信息;批检与单检共用) */
  function searchPlatformFor(g: string): { platform: string; type: string } {
    const row = searchResults.value.find((r) => String(r.gse ?? '').toUpperCase() === g)
    return { platform: String(row?.platform ?? ''), type: String(row?.type ?? '') }
  }

  async function probeBatch(force = false): Promise<void> {
    pruneBadges()
    const rows = searchResults.value
      .filter((r) => selectedRows.value.has(rowGse(r)))
      // 已有新鲜结论(TTL 内)的行跳过,重复打网络无增量;过期/err 行自动重测;
      // Shift+点击 = 强制重测(忽略缓存,QA 2026-11-09),本地数据集状态变化
      // 另有失效钩子(runJob/forceJob onDone)精准作废对应条目(总览不触发)
      .filter((r) => force || !freshBadge(rowGse(r)))
      .map((r) => ({ gse: rowGse(r), ...searchPlatformFor(rowGse(r)) }))
    if (!rows.length) { flash(t('web.ml.probe_sel_cached')); return }
    probeBatchBusy.value = true
    try {
      const res = await api.mlProbeBatch(rows)
      for (const item of res.results) putBadge(item)
      flash(t('web.ml.probe_batch_done', { n: res.results.length }))
    } catch (e) {
      flash(t('web.ml.probe_fail', { msg: e instanceof Error ? e.message : String(e) }))
    } finally {
      probeBatchBusy.value = false
    }
  }

  return {
    selectedRows, probeBatchBusy, selectedCount, allSelected,
    rowGse, toggleRow, toggleAll, symbolVerdict, dotClass,
    searchPlatformFor, probeBatch
  }
}
