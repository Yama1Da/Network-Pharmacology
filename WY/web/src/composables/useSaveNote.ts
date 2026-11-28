import { ref, watch, type WatchSource } from 'vue'
import { useI18n } from 'vue-i18n'

/** 保存结果内联提示(SettingsView 四栏共用):页面长,顶部 notice 在内容区外
 * 看不见,随保存按钮就地显示;切换板块即清空,避免上一板块的结果串场。
 * transfer 的导出/导入流程仍走顶部 notice(useNotice),不经此路。
 * resetOn 传入"当前板块"响应源,切换即清空。在 setup 上下文中调用。 */
export function useSaveNote(resetOn?: WatchSource) {
  const { t } = useI18n()
  const saveNote = ref<{ ok: boolean; text: string } | null>(null)
  let saveNoteTimer = 0

  function noteFlash(ok: boolean, text: string): void {
    saveNote.value = { ok, text }
    window.clearTimeout(saveNoteTimer)
    saveNoteTimer = window.setTimeout(() => (saveNote.value = null), 4000)
  }

  /** 保存后的统一提示:基础文案 + 整合缓存清理说明(随保存按钮就地显示) */
  function flashSaved(purged: number): void {
    noteFlash(true, purged > 0
      ? `${t('web.settings.saved')} ${t('web.settings.purged', { n: purged })}`
      : t('web.settings.saved'))
  }

  if (resetOn) watch(resetOn, () => (saveNote.value = null))

  return { saveNote, noteFlash, flashSaved }
}
