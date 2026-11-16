import { computed, ref, watch } from 'vue'

/** 列表分页:filtered 惰性传入(调用方的过滤 computed 通常依赖各自的搜索词,
 * 二者存在定义顺序上的相互引用);goto 自动夹紧边界。
 * nPages 收缩(删除/后台刷新把条目刷少)时自动把 page 拉回有效范围:
 * 不收紧的话停在超出页会看到空列表 + "5 / 4" 页码(2026-11 审计 M1)。 */
export function usePager<T>(filtered: () => T[], pageSize = 10) {
  const page = ref(1)
  const paged = computed(() =>
    filtered().slice((page.value - 1) * pageSize, page.value * pageSize))
  const nPages = computed(() => Math.max(1, Math.ceil(filtered().length / pageSize)))
  watch(nPages, () => {
    if (page.value > nPages.value) page.value = nPages.value
  })
  function goto(p: number): void {
    page.value = Math.min(Math.max(1, p), nPages.value)
  }
  return { page, paged, nPages, goto }
}
