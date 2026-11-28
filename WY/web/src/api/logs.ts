/** 日志(阶段 1):日志文件列表/内容读取与级别调整,端点与 api/logs.py 一一对应。 */

import { json, post } from './http'

export interface LogFile {
  name: string
  modified: string
  size: number
  is_current: boolean
}

interface LogContent {
  name: string
  content: string
  truncated: boolean
  size: number
}

export const logsApi = {
  logsFiles: (): Promise<{ files: LogFile[] }> =>
    fetch('/api/logs/files').then((r) => json(r)),

  logsContent: (name: string, maxBytes = 2_000_000): Promise<LogContent> =>
    fetch(`/api/logs/content?name=${encodeURIComponent(name)}&max_bytes=${maxBytes}`).then(
      (r) => json(r)
    ),

  logsLevel: (): Promise<{ level: string }> =>
    fetch('/api/logs/level').then((r) => json(r)),

  setLogsLevel: (level: 'INFO' | 'DEBUG'): Promise<{ level: string }> =>
    post('/api/logs/level', { level })
}
