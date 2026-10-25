import { createI18n } from 'vue-i18n'
import en from '../locales/en.json'
import zh from '../locales/zh.json'

/** vue-i18n:消息源为 scripts/export_i18n_json.py 导出的词典(勿手改,重跑即刷新)。
 * 语言偏好存 localStorage('npap-lang');默认中文。 */
export const i18n = createI18n({
  legacy: false,
  locale: localStorage.getItem('npap-lang') || 'zh',
  fallbackLocale: 'zh',
  messages: { zh, en }
})

export function setLang(lang: 'zh' | 'en'): void {
  i18n.global.locale.value = lang
  localStorage.setItem('npap-lang', lang)
}
