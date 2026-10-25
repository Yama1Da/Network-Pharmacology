import { onBeforeUnmount, ref } from 'vue'

/** 两段式确认:首次点击进入武装态(timeout 毫秒内未再点自动撤防),
 * 武装态下再点执行 confirm 并撤防。模板以 armed 驱动按钮的 arm 样式与文案切换。 */
export function useConfirmArm(timeout = 4000) {
  const armed = ref(false)
  let timer = 0
  function tap(confirm: () => void): void {
    if (!armed.value) {
      armed.value = true
      window.clearTimeout(timer)
      timer = window.setTimeout(() => (armed.value = false), timeout)
      return
    }
    armed.value = false
    window.clearTimeout(timer)
    confirm()
  }
  onBeforeUnmount(() => window.clearTimeout(timer))
  return { armed, tap }
}
