import { reactive } from 'vue'

/** 跨视图轻量状态:主页检索完成后向分析页预填药物/疾病(单页应用内共享,无持久化)。 */
export const prefill = reactive<{ drug: string; disease: string; stamp: number }>({
  drug: '',
  disease: '',
  stamp: 0
})

export function setPrefill(drug: string, disease: string): void {
  prefill.drug = drug
  prefill.disease = disease
  prefill.stamp = Date.now()
}

export function takePrefill(): { drug: string; disease: string } | null {
  if (!prefill.stamp) return null
  const out = { drug: prefill.drug, disease: prefill.disease }
  prefill.stamp = 0
  return out
}

/** 跨视图轻量状态:协同分析结果页「复现此运行」→ 协同预测页回填表单与参数。 */
export interface SynergyReplay {
  drug_a: string
  drug_b: string
  disease: string
  top_n: number | null
  n_perm: number | null
}

export const synergyReplay = reactive<{ payload: SynergyReplay | null; stamp: number }>({
  payload: null,
  stamp: 0
})

export function setSynergyReplay(payload: SynergyReplay): void {
  synergyReplay.payload = payload
  synergyReplay.stamp = Date.now()
}

export function takeSynergyReplay(): SynergyReplay | null {
  if (!synergyReplay.stamp) return null
  const out = synergyReplay.payload
  synergyReplay.stamp = 0
  return out
}
