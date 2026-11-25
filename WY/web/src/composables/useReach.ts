/**
 * 逐靶点可达性标注(W7,2026-11-25):主页三表与分析页 hub 表共用的懒取 +
 * 过滤逻辑。标注列默认开(数据到了就渲染,取数失败静默降级=列不出现,不影响
 * 一期结果);"过滤可关、默认只标注"(方案 §答辩FAQ)——过滤开关默认"全部"
 * 保一期行为。"暴露未知"独立显示,不与"不可达"合并(四档词表复用
 * web.synergy.reach_tier_*,单一真相源)。
 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type TargetsReachability } from '../api'

export type ReachFilter = 'all' | 'partial_up' | 'reachable' | 'judged'

export function useReachAnnotation() {
  const { t } = useI18n()
  const reach = ref<TargetsReachability | null>(null)
  const reachFilter = ref<ReachFilter>('all')

  /** 懒取指定药物的逐靶点档位;不阻塞调用方主结果渲染 */
  async function load(drug: string): Promise<void> {
    reach.value = null
    reachFilter.value = 'all'
    const name = (drug || '').trim()
    if (!name) return
    try {
      reach.value = await api.targetsReachability(name)
    } catch {
      /* 静默降级:标注列不出现 */
    }
  }

  function reset(): void {
    reach.value = null
    reachFilter.value = 'all'
  }

  const hasReach = computed(() =>
    !!reach.value?.reachability && Object.keys(reach.value.genes).length > 0)

  /** 基因 → 四档档位(无标注数据时归 unknown,列不渲染故不影响展示) */
  function tier(gene: unknown): string {
    return reach.value?.genes[String(gene ?? '')]?.tier ?? 'unknown'
  }

  /** unknown 档按真因拆词(2026-11-25 负责人定稿):分子缺失(no_exposure)=
   * "暴露未知",分母缺失(no_activity)= "活性未知"——metformin 暴露已知而
   * 靶点无活性记录的行,称"暴露未知"是误称;两词仍与"不可达"严格区分 */
  function word(gene: unknown): string {
    const g = reach.value?.genes[String(gene ?? '')]
    if (!g) return t('web.synergy.reach_tier_unknown')
    if (g.tier === 'unknown' && g.unknown_reason === 'no_activity')
      return t('web.reach.tier_no_activity')
    return t(`web.synergy.reach_tier_${g.tier}`)
  }

  /** 悬浮注解:R 值 + 档位词(R 为空只给档位词;unknown 行附真因词) */
  function title(gene: unknown): string {
    const w = word(gene)
    const g = reach.value?.genes[String(gene ?? '')]
    if (g?.r_value == null) return w
    return `R=${Number(g.r_value).toPrecision(3)} · ${w}`
  }

  /** 过滤语义:partial_up=可达+部分;reachable=仅可达;judged=已判档(排除未知) */
  function matchesFilter(tierName: string, f: ReachFilter): boolean {
    if (f === 'all') return true
    if (f === 'reachable') return tierName === 'reachable'
    if (f === 'partial_up') return tierName === 'reachable' || tierName === 'partial'
    return tierName !== 'unknown'
  }

  return { reach, reachFilter, load, reset, hasReach, tier, word, title, matchesFilter }
}
