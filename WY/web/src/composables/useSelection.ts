import { computed, ref } from 'vue'

/** 勾选集合:单项切换 / 全选(allNames 口径与既有页面一致,按全量而非过滤结果)/
 * 数据重载后 retain() 剔除已消失的名字,避免悬空勾选。 */
export function useSelection(allNames: () => Iterable<string>) {
  const selected = ref<Set<string>>(new Set())
  function toggle(name: string): void {
    const next = new Set(selected.value)
    if (next.has(name)) next.delete(name)
    else next.add(name)
    selected.value = next
  }
  const allSelected = computed(() => {
    const all = [...allNames()]
    return all.length > 0 && all.every((n) => selected.value.has(n))
  })
  function toggleAll(): void {
    selected.value = allSelected.value ? new Set() : new Set(allNames())
  }
  function retain(names: Iterable<string>): void {
    const keep = new Set(names)
    selected.value = new Set([...selected.value].filter((n) => keep.has(n)))
  }
  return { selected, toggle, allSelected, toggleAll, retain }
}
