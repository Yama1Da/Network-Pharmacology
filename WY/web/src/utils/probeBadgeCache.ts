/** ML 页 GSE 体检徽章的会话级缓存(QA 2026-11-09)。
 *
 * 徽章判定(matrix_empty / counts / 注释)是 GEO 归档的不可变属性,缓存语义安全;
 * 唯一会翻转它的是本地盘上状态(数据集下载/强制重下完成)——由视图在对应 job 的
 * onDone 里调 dropBadge() 精准失效;TTL 1h 兜底,过期自动重测。
 * 纯浏览器内存,不落盘,刷新页面即清空;容量上限防多疾病检索累积。
 * err(网络层瞬时失败)不入缓存:下次批检自然重试,免得把瞬时故障钉死一个 TTL。
 */
import type { MlProbeBatchItem } from '../api'

const TTL_MS = 60 * 60 * 1000
const MAX_ENTRIES = 300

interface Entry {
  item: MlProbeBatchItem
  at: number
}

const cache = new Map<string, Entry>()

/** 未过期的徽章结论;过期即删并按不存在处理(灰点,批检重测)。 */
export function freshBadge(gse: string): MlProbeBatchItem | undefined {
  const e = cache.get(gse)
  if (!e) return undefined
  if (Date.now() - e.at > TTL_MS) {
    cache.delete(gse)
    return undefined
  }
  return e.item
}

export function putBadge(item: MlProbeBatchItem): void {
  if (item.status === 'err') return
  cache.set(item.gse, { item, at: Date.now() })
  while (cache.size > MAX_ENTRIES) {
    const oldest = cache.keys().next().value
    if (oldest === undefined) break
    cache.delete(oldest)
  }
}

/** 失效钩子:数据集下载/强制重下完成后由视图调用,该 GSE 徽章立刻回灰点重测。 */
export function dropBadge(gse: string): void {
  cache.delete(gse)
}

/** 剪掉过期条目(检索/批检入口调用即可;过期语义另由 freshBadge 逐条兜底)。 */
export function pruneBadges(): void {
  const now = Date.now()
  for (const [g, e] of cache) {
    if (now - e.at > TTL_MS) cache.delete(g)
  }
}
