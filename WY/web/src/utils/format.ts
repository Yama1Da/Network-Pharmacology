/** 共享格式化:文件大小(Logs/Sources/Collection/Ml 四页同款)与数值列(Analysis/Docking 同款)。 */
export const fmtSize = (n: number): string =>
  n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(1)} KB` : `${(n / 1048576).toFixed(1)} MB`

export const fmt = (v: unknown, digits = 3): string =>
  typeof v === 'number' ? v.toFixed(digits) : (v == null ? '—' : String(v))
