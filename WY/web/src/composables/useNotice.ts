import { onBeforeUnmount, ref } from 'vue'

/** 轻提示:flash() 写入一条 notice,timeout 毫秒后自动消失;组件卸载时清理定时器。
 * 各视图模板直接消费 notice(配 Transition name="note")。 */
export function useNotice(timeout = 3200) {
  const notice = ref('')
  let timer = 0
  function flash(msg: string): void {
    notice.value = msg
    window.clearTimeout(timer)
    timer = window.setTimeout(() => (notice.value = ''), timeout)
  }
  onBeforeUnmount(() => window.clearTimeout(timer))
  return { notice, flash }
}
