import { computed, ref } from 'vue'

/** 列表分页:filtered 惰性传入(调用方的过滤 computed 通常依赖各自的搜索词,
 * 二者存在定义顺序上的相互引用);goto 自动夹紧边界。 */
export function usePager<T>(filtered: () => T[], pageSize = 10) {
  const page = ref(1)
  const paged = computed(() =>
    filtered().slice((page.value - 1) * pageSize, page.value * pageSize))
  const nPages = computed(() => Math.max(1, Math.ceil(filtered().length / pageSize)))
  function goto(p: number): void {
    page.value = Math.min(Math.max(1, p), nPages.value)
  }
  return { page, paged, nPages, goto }
}
