/** 界面等比缩放(设置页·通用设置的"字体大小"):对文档根元素做 CSS zoom,
 * 文字按缩放后尺寸重新光栅化(不发虚),布局随缩放自动重排。
 * 双层持久化:localStorage 供启动即时恢复(避免先 100% 再跳变),
 * 后端 ui_settings(UI_FONT_SCALE)为权威值,App 启动后异步校准。 */

export const FONT_SCALE_MIN = 0.8
export const FONT_SCALE_MAX = 1.5
export const FONT_SCALE_STEP = 0.1

const LS_KEY = 'npap_font_scale'

export function clampScale(v: unknown): number {
  const n = Number(v)
  if (!Number.isFinite(n)) return 1
  return Math.min(FONT_SCALE_MAX, Math.max(FONT_SCALE_MIN, Math.round(n * 100) / 100))
}

/** 吸附到步长栅格(旧档位 0.85/1.25 等不在 10% 栅格上的值就近取整),并清掉浮点尾数 */
export function snapScale(v: unknown): number {
  const c = clampScale(v)
  return Math.round((Math.round(c / FONT_SCALE_STEP) * FONT_SCALE_STEP) * 100) / 100
}

export function storedScale(): number {
  try {
    const raw = window.localStorage.getItem(LS_KEY)
    return raw == null ? 1 : clampScale(raw)
  } catch {
    return 1
  }
}

export function persistScale(v: number): void {
  try {
    window.localStorage.setItem(LS_KEY, String(clampScale(v)))
  } catch {
    /* 隐私模式等 localStorage 不可用时放弃本地镜像,仅后端持久化 */
  }
}

type ZoomableStyle = CSSStyleDeclaration & { zoom?: string }

export function applyScale(v: number): void {
  const s = clampScale(v)
  const style = document.documentElement.style as ZoomableStyle
  style.zoom = s === 1 ? '' : String(s)  // 空串 = 移除内联 zoom,回到默认 100%
}
