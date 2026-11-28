/** 跨域共享类型:api/ 各域模块与组装层 api.ts 共用,只放类型不放运行时逻辑。 */

export type Row = Record<string, unknown>

export interface JobState {
  job_id: string
  status: 'running' | 'done' | 'error'
  description: string
  progress: { fraction: number | null; text: string; [k: string]: unknown }
  result: Record<string, unknown> | null
  error: string | null
  finished?: number | null   // 完成时刻(秒);回填窗口与任务中心倒计时以此为准
  trace?: string
}

/** 在途任务状态(/api/{ml,ppi,synergy}/job 的公共形态):重进页面恢复进度用。 */
export interface InFlightJob {
  running: boolean
  job_id?: string
  status?: string
  progress?: JobState['progress']
}
