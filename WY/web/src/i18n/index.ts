import { createI18n } from 'vue-i18n'
import en from '../locales/en.json'
import zh from '../locales/zh.json'

export type Lang = 'zh' | 'en'

const LS_KEY = 'npap-lang'

export function normalizeLang(v: unknown): Lang {
  return v === 'en' ? 'en' : 'zh'
}

export function storedLang(): Lang {
  try {
    const raw = window.localStorage.getItem(LS_KEY)
    return raw == null ? 'zh' : normalizeLang(raw)
  } catch {
    return 'zh'
  }
}

export function persistLang(lang: Lang): void {
  try {
    window.localStorage.setItem(LS_KEY, lang)
  } catch {
    /* localStorage 不可用时放弃本地镜像,仅后端持久化 */
  }
}

/** vue-i18n:消息源为 scripts/export_i18n_json.py 导出的词典(勿手改,重跑即刷新)。
 * 语言偏好双层持久化:localStorage 供启动即时恢复(避免闪切换),后端 UI_LANG 为权威值,
 * App 启动后异步校准,并随 ui_settings.json 参与数据迁移(与颜色主题同机制)。 */
export const i18n = createI18n({
  legacy: false,
  locale: storedLang(),
  fallbackLocale: 'zh',
  messages: { zh, en }
})

export function setLang(lang: Lang): void {
  i18n.global.locale.value = normalizeLang(lang)
  persistLang(lang)
}
