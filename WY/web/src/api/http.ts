/** fetch 解包与请求原语:api/ 各域模块与组装层 api.ts 共用。 */

/** 统一响应解包:失败时优先取 FastAPI 的 detail 文本(设置校验、来源校验等都
 * 写成了人话,只在状态码上断掉会让用户看到裸 "HTTP 400"),取不到再退回状态码
 * (2026-11 审计 M5:原 json() 丢弃 detail,全端点 4xx 只剩 "HTTP 400")。 */
export async function json<T>(r: Response): Promise<T> {
  if (!r.ok) {
    let detail = ''
    try {
      detail = ((await r.json()) as { detail?: string })?.detail ?? ''
    } catch {
      /* 响应体不是 JSON 时退回状态码 */
    }
    throw new Error(detail || `HTTP ${r.status}`)
  }
  return (await r.json()) as T
}

/** 兼容别名:detail 提取已并入 json(),所有端点统一受益。 */
export const jsonWithDetail = json

export async function post<T>(url: string, body: unknown): Promise<T> {
  return fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body)
  }).then((r) => json<T>(r))
}

/** 从 Content-Disposition 取服务端文件名:优先 RFC 5987 的 filename*(中文药名场景),
 * 退回普通 filename="…";都缺失返回 null,由调用方兜底命名。 */
export function attachmentName(disposition: string | null): string | null {
  if (!disposition) return null
  const star = /filename\*\s*=\s*(?:UTF-8|utf-8)''([^;]+)/.exec(disposition)
  if (star) {
    try {
      return decodeURIComponent(star[1].trim())
    } catch {
      /* 非法百分号编码时退回普通 filename */
    }
  }
  const plain = /filename\s*=\s*"([^"]*)"/.exec(disposition) ??
                /filename\s*=\s*([^;]+)/.exec(disposition)
  return plain ? plain[1].trim() : null
}
