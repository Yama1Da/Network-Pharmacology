/** NPAP API 的类型化访问层:与 WY/api.py 的端点一一对应。 */

export interface LogFile {
  name: string
  modified: string
  size: number
  is_current: boolean
}

export interface LogContent {
  name: string
  content: string
  truncated: boolean
  size: number
}

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) throw new Error(`HTTP ${r.status}`)
  return (await r.json()) as T
}

export const api = {
  files: (): Promise<{ files: LogFile[] }> =>
    fetch('/api/logs/files').then((r) => json(r)),

  content: (name: string, maxBytes = 2_000_000): Promise<LogContent> =>
    fetch(`/api/logs/content?name=${encodeURIComponent(name)}&max_bytes=${maxBytes}`).then(
      (r) => json(r)
    ),

  level: (): Promise<{ level: string }> =>
    fetch('/api/logs/level').then((r) => json(r)),

  setLevel: (level: 'INFO' | 'DEBUG'): Promise<{ level: string }> =>
    fetch('/api/logs/level', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ level })
    }).then((r) => json(r))
}
