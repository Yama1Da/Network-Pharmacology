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
