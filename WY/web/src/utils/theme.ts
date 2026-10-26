/** 颜色主题(设置页·通用设置的"颜色变换"):浅色/深色(样品「墨黑」)。
 * 切换 = html[data-theme],调色板见 screen.css 顶部变量表。
 * 双层持久化:localStorage 供启动即时恢复(避免闪浅色),后端 UI_THEME 为权威值,
 * App 启动后异步校准(与字体大小设置同机制)。theme 为响应式 ref,供 venn 图等随动。 */
import { ref } from 'vue'

export type Theme = 'light' | 'dark'

export const theme = ref<Theme>('light')

const LS_KEY = 'npap_theme'

export function normalizeTheme(v: unknown): Theme {
  return v === 'dark' ? 'dark' : 'light'
}

export function storedTheme(): Theme {
  try {
    const raw = window.localStorage.getItem(LS_KEY)
    return raw == null ? 'light' : normalizeTheme(raw)
  } catch {
    return 'light'
  }
}

export function persistTheme(t: Theme): void {
  try {
    window.localStorage.setItem(LS_KEY, t)
  } catch {
    /* localStorage 不可用时放弃本地镜像,仅后端持久化 */
  }
}

export function applyTheme(t: Theme): void {
  theme.value = normalizeTheme(t)
  if (theme.value === 'dark') document.documentElement.dataset.theme = 'dark'
  else delete document.documentElement.dataset.theme
}
