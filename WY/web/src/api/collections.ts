/** 文件集合(阶段 2)与来源刷新(管理页):downloads 各类产物的列举/打包/删除
 * 与 drug/disease 来源刷新;kind 枚举与后端集合注册表一致。 */

import { json, post } from './http'

export type CollectionKind =
  | 'drug_sources' | 'disease_sources'
  | 'ppi' | 'enrichment' | 'stars' | 'downloads'

export interface CollectionItem {
  name: string
  modified: string
  size: number
  /* 管理页(drug/disease_sources)附加 */
  is_merged?: boolean
  merged_status?: 'current' | 'legacy' | null
  source?: string
  /* 富集页附加 */
  display?: string
  side?: string
  gsea?: boolean
}

interface CollectionListing {
  kind: CollectionKind
  mode: 'files' | 'dirs'
  items: CollectionItem[]
}

interface DeleteResult {
  deleted: number
  failed: string[]
}

export const collectionsApi = {
  collection: (kind: CollectionKind): Promise<CollectionListing> =>
    fetch(`/api/collections/${kind}`).then((r) => json(r)),

  zip: async (kind: CollectionKind, names: string[]): Promise<Blob> => {
    const r = await fetch(`/api/collections/${kind}/zip`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ names })
    })
    if (!r.ok) throw new Error(`HTTP ${r.status}`)
    return r.blob()
  },

  remove: (kind: CollectionKind, names: string[]): Promise<DeleteResult> =>
    post(`/api/collections/${kind}/delete`, { names }),

  /* ---- 来源刷新(管理页) ---- */
  refreshSources: (
    kind: 'drug' | 'disease',
    filenames: string[] | null,
    hours?: number
  ): Promise<{ job_id: string }> =>
    post(`/api/sources/${kind}/refresh`, { filenames, hours })
}
