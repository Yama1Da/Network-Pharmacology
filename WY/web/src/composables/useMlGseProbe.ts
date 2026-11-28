import { computed, ref, type Ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type MlProbeResult } from '../api'
import { freshBadgeWithAge } from '../utils/probeBadgeCache'

/** GSE 单检与可用性锁定(MlView ① 区)。GSE 确定后预检(每 GSE 一次;
 * 徽章命中则豁免,2026-11-10 QA 拍板):批检与单检同走后端 _probe_one,
 * 1h TTL 内绿/红结论相同——命中即跳过单检(省一次 Range 请求与缓存命中时的
 * 全矩阵解析),红点直接锁定;灰(unknown)、err(不缓存,重试有真实价值)、
 * 无徽章(直填/超时)照旧单检。unusable 是②③两区的禁用总闸(徽章直锁 +
 * assess 判定)。本地矩阵三输入(expr/labels/case)也在此——与 GSE 构成
 * 互斥输入源。在 setup 上下文中调用。 */
export function useMlGseProbe(
  gse: Ref<string>,
  searchPlatformFor: (g: string) => { platform: string; type: string }
) {
  const { t } = useI18n()
  const exprPath = ref('')
  const labelsPath = ref('')
  const localCase = ref('')
  const probe = ref<MlProbeResult | null>(null)
  const probeBusy = ref(false)
  let probedGse = ''
  const badgeApplied = ref<{ status: 'ok' | 'unusable'; ageMin: number } | null>(null)

  async function probeGse(): Promise<void> {
    const g = gse.value.trim().toUpperCase()
    badgeApplied.value = null
    if (!g) return
    // 徽章查询必须在"每 GSE 一次"去重之前:批检徽章常在单检已跑完后才落地
    // (黄点期间点过该 GSE),再次 change 同一 GSE 仍要查——否则占位说明永远
    // 不会出现,单检的旧结果(可能是 err)一直占着展示位
    const b = freshBadgeWithAge(g)
    if (b && (b.item.status === 'ok' || b.item.status === 'unusable')) {
      probedGse = g
      badgeApplied.value = { status: b.item.status, ageMin: b.ageMin }
      probe.value = null
      probeBusy.value = false
      return
    }
    if (probedGse === g) return
    probedGse = g
    probeBusy.value = true
    try {
      probe.value = await api.mlProbe(g, searchPlatformFor(g))
    } catch (e) {
      probe.value = { err: e instanceof Error ? e.message : String(e) }
    } finally {
      probeBusy.value = false
    }
  }

  function onGseChange(): void {
    probedGse = ''
    probe.value = null
    void probeGse()
  }

  /** 重置单检状态(检索新词 / 强制重下落盘后调用) */
  function resetProbe(): void {
    probedGse = ''
    probe.value = null
    badgeApplied.value = null
  }

  const useLocal = computed(() => !!exprPath.value.trim() && !!labelsPath.value.trim())
  const localPartial = computed(() => !!exprPath.value.trim() !== !!labelsPath.value.trim())
  const hasGse = computed(() => !!gse.value.trim())
  const sourceConflict = computed(() => useLocal.value && hasGse.value)

  /* 预检展示与可用性锁定 */
  const probeR = computed(() => probe.value?.r ?? null)
  const probeErr = computed(() => probe.value?.err ?? '')
  const probeAssess = computed(() => probe.value?.assess ?? null)
  const unusable = computed(() => {
    // 徽章直锁:批检判 unusable 的红点数据集跳过单检,选中即锁(2026-11-10 QA;
    // 批检返回体无 annot_missing 字段,不能靠合成 assess 判定,须按 status 直锁)
    if (hasGse.value && badgeApplied.value?.status === 'unusable') return true
    const p = probe.value
    const a = probeAssess.value
    return !!(hasGse.value && a && p?.gse === gse.value.trim().toUpperCase() &&
      ((a.matrix_empty && !a.counts_available) ||
       (!a.symbol_level && a.annot === 'missing')))
  })

  const fmtEta = (s: number): string =>
    s < 90 ? `${Math.max(1, Math.round(s))} ${t('web.ml.unit_sec')}`
      : s < 3600 ? `${Math.max(1, Math.round(s / 60))} ${t('web.ml.unit_min')}`
        : `${Math.floor(s / 3600)} ${t('web.ml.unit_hour')} ${Math.round(s % 3600 / 60)} ${t('web.ml.unit_min_short')}`

  return {
    exprPath, labelsPath, localCase,
    probe, probeBusy, badgeApplied,
    probeGse, onGseChange, resetProbe,
    useLocal, localPartial, hasGse, sourceConflict,
    probeR, probeErr, probeAssess, unusable, fmtEta
  }
}
