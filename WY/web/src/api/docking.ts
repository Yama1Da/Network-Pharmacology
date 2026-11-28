/** 分子对接(阶段 3b):sid 会话驱动的多步流程(source → 结构/受体 → box → run → save),
 * 端点与 api/docking.py 一一对应;载荷类型内联在方法签名上(无独立具名类型)。 */

import { json, jsonWithDetail, post } from './http'
import type { Row } from './types'

export const dockingApi = {
  dockingSession: (): Promise<{ sid: string }> =>
    post('/api/docking/session', {}),

  dockingSources: (): Promise<{
    ppi_dirs: { name: string; display: string }[]
    enr_files: { name: string; display: string }[]
    manual_dirs: { name: string; n_genes: number; mtime: string; genes: string[] }[]
  }> => fetch('/api/docking/sources').then((r) => json(r)),

  dockingSource: (
    sid: string,
    body: { mode: 'ppi' | 'enr' | 'manual' | 'manual-file'; name?: string; manual_genes?: string; manual_file?: string }
  ): Promise<{ label: string; drug: string; disease: string; hub_count: number;
               enr_count: number; candidates: string[] }> =>
    post(`/api/docking/${sid}/source`,
         { mode: body.mode, name: body.name ?? '', manual_genes: body.manual_genes ?? '',
           manual_file: body.manual_file ?? '' }),

  dockingCheckStructures: (
    sid: string,
    body: { use_struct: boolean; allow_af: boolean }
  ): Promise<{ job_id: string | null; direct?: Record<string, unknown> }> =>
    post(`/api/docking/${sid}/check-structures`, body),

  dockingLigand: (
    sid: string,
    body: { source: 'pubchem' | 'smiles'; query?: string; smiles?: string }
  ): Promise<{ ok: boolean; error?: string; cid?: number | null; lig_pdbqt?: string }> =>
    post(`/api/docking/${sid}/ligand`, body),

  dockingReceptorCandidates: (
    sid: string,
    body: { gene: string; allow_af: boolean }
  ): Promise<{ error?: string; rows: { pdb_id: string; resolution: string; method: string;
             title: string; source: string; flag: string; note: string }[]; af_enabled?: boolean }> =>
    post(`/api/docking/${sid}/receptor/candidates`, body),

  dockingReceptorPrepare: (
    sid: string,
    body: { gene: string; pdb_id: string }
  ): Promise<{ job_id: string }> =>
    post(`/api/docking/${sid}/receptor/prepare`, body),

  // 自带受体上传:原始请求体直传(仿 transferImport,不引 multipart),后端
  // detail(扩展名/大小/体检硬伤)经 jsonWithDetail 直达 UI
  dockingReceptorUpload: (sid: string, gene: string, file: File): Promise<{ job_id: string }> =>
    fetch(
      `/api/docking/${sid}/receptor/upload?gene=${encodeURIComponent(gene)}&name=${encodeURIComponent(file.name)}`,
      { method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: file }
    ).then((r) => jsonWithDetail<{ job_id: string }>(r)),

  dockingBox: (
    sid: string,
    body: { infer: boolean; center?: number[]; size?: number[] }
  ): Promise<{ box: { center: number[]; size: number[] }; clamped: boolean;
               inferred_ligands: string[] | null; residues: string[] }> =>
    post(`/api/docking/${sid}/box`, body),

  dockingRun: (sid: string): Promise<{ job_id: string }> =>
    post(`/api/docking/${sid}/run`, {}),

  dockingSave: (
    sid: string,
    body: { pose_mode: number; conclusion: string; drug: string; disease: string }
  ): Promise<{ result_dir: string }> =>
    post(`/api/docking/${sid}/save`, body),

  dockingStar: (sid: string, on: boolean): Promise<{ starred: boolean }> =>
    post(`/api/docking/${sid}/star`, { on }),

  dockingInteractions: (sid: string, pose: number): Promise<{ rows: Row[] }> =>
    fetch(`/api/docking/${sid}/interactions?pose=${pose}`).then((r) => json(r)),

  dockingView3dUrl: (sid: string, pose: number, cutoff: number, focus: boolean): string =>
    `/api/docking/${sid}/view3d?pose=${pose}&cutoff=${cutoff}&focus=${focus ? 1 : 0}`,

  dockingSave3d: (
    sid: string,
    body: { pose: number; caption: string; cutoff: number }
  ): Promise<{ path: string; artifact_rel: string; offline_ok: boolean }> =>
    post(`/api/docking/${sid}/save3d`, body)
}
