<script setup lang="ts">
/** 主页(App.py 职责的 TS 化):四库检索落地页,1:1 复刻 Front-end-screen-options 极简主义样本。
 * 箭头点击 → /api/targets/retrieve job(四库并发检索,写入来源与 merged 缓存),
 * 完成后展示两侧靶点数并可预填跳转分析页。 */
import { onBeforeUnmount, onMounted, ref, type Ref } from 'vue'
import { api } from '../api'
import { setPrefill } from '../store'
import { useJob } from '../composables/useJob'

const drug = ref('metformin')
const disease = ref('Type 2 diabetes mellitus')

const job = useJob((result) => {
  done.value = {
    drug: Number(result.drug_count ?? 0),
    disease: Number(result.disease_count ?? 0),
    common_count: Number(result.common_count ?? 0),
    common_genes: (result.common_genes ?? []) as string[],
    failed: (result.failed ?? []) as string[]
  }
})

const done = ref<{
  drug: number
  disease: number
  common_count: number
  common_genes: string[]
  failed: string[]
} | null>(null)

function retrieve(): void {
  if (!drug.value.trim() || !disease.value.trim()) {
    return
  }
  done.value = null
  void job.run(api.retrieveTargets(drug.value.trim(), disease.value.trim()),
    { doneText: '四库检索完成' })
}

function goAnalysis(): void {
  setPrefill(drug.value.trim(), disease.value.trim())
  location.hash = '#/analysis'
}

/* 页脚统计数字滚动(样本同款:装饰性品牌数字) */
interface CountStat {
  v: Ref<number>
  t: number
}
const stats: CountStat[] = [
  { v: ref(0), t: 12 },
  { v: ref(0), t: 8932 },
  { v: ref(0), t: 31486 },
  { v: ref(0), t: 402 }
]
const fmt = (n: number): string => n.toLocaleString('en-US')

let raf = 0
let started = false
function startCount(): void {
  if (started) return
  started = true
  const t0 = performance.now()
  const dur = 1600
  const step = (t: number): void => {
    const p = Math.min(1, (t - t0) / dur)
    const e = 1 - Math.pow(1 - p, 3)
    for (const s of stats) s.v.value = Math.round(s.t * e)
    if (p < 1) raf = window.requestAnimationFrame(step)
  }
  raf = window.requestAnimationFrame(step)
}
onMounted(() => window.setTimeout(startCount, 1750))  // 与 App 入场幕布对齐
onBeforeUnmount(() => window.cancelAnimationFrame(raf))
</script>

<template>
  <p class="kicker s1">NETWORK&nbsp;PHARMACOLOGY</p>
  <h1 class="title s2">网络药理学分析平台</h1>
  <p class="sub s3">成分&nbsp;—&nbsp;靶点&nbsp;—&nbsp;通路&nbsp;—&nbsp;验证</p>

  <div class="home-form s4">
    <label class="field">
      <input v-model="drug" spellcheck="false" placeholder="药物成分" />
      <i class="line"></i>
    </label>
    <span class="cross">×</span>
    <label class="field">
      <input v-model="disease" spellcheck="false" placeholder="疾病表型" />
      <i class="line"></i>
    </label>
    <button class="go" :disabled="job.running.value" title="四库检索" aria-label="四库检索"
            @click="retrieve">
      <svg v-if="!job.running.value" viewBox="0 0 24 14" width="22" height="13">
        <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
              stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
      </svg>
      <span v-else class="spin"></span>
    </button>
  </div>

  <Transition name="note">
    <div v-if="job.running.value" class="retrieve">
      <div class="progress">
        <div class="progress-track">
          <div class="progress-fill"
               :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
        </div>
        <p class="progress-text">{{ job.text.value || '四库检索中…' }}</p>
      </div>
    </div>
  </Transition>

  <p v-if="job.error.value" class="retrieve failed-msg">
    检索失败:{{ job.error.value }} —— 若刚更新过程序,请重启 npap serve 并强制刷新(Ctrl+F5)后重试。
  </p>

  <Transition name="note">
    <div v-if="done" class="retrieve done">
      <p class="done-line">
        药物 <b>{{ done.drug }}</b> 靶点 · 疾病 <b>{{ done.disease }}</b> 靶点 ·
        交集 <b>{{ done.common_count }}</b>
        <span v-if="done.failed.length" class="failed">(失败来源:{{ done.failed.join(', ') }})</span>
      </p>
      <img v-if="job.jobId.value" class="venn" :src="`/api/targets/venn?job_id=${job.jobId.value}`"
           alt="药物-疾病靶点 Venn 图" />
      <p v-if="done.common_genes.length" class="hint-line">
        交集靶点:{{ done.common_genes.slice(0, 24).join(', ') }}
        {{ done.common_genes.length > 24 ? `等 ${done.common_genes.length} 个` : '' }}
      </p>
      <button class="go enter" title="进入分析" @click="goAnalysis">
        <svg viewBox="0 0 24 14" width="22" height="13">
          <path d="M1 7h20M15.5 1.5 21.5 7l-6 5.5" fill="none" stroke="currentColor"
                stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" />
        </svg>
      </button>
      <p class="hint-line">来源与整合缓存已就绪,分析页将直接命中。</p>
    </div>
  </Transition>

  <div class="foot-stats s5">
    <span><b>{{ fmt(stats[0].v.value) }}</b> 数据库 ·
      <b>{{ fmt(stats[1].v.value) }}</b> 化合物 ·
      <b>{{ fmt(stats[2].v.value) }}</b> 靶点 ·
      <b>{{ fmt(stats[3].v.value) }}</b> 通路</span>
  </div>
</template>

<style scoped>
.home-form {
  margin-top: 56px;
  display: flex;
  align-items: flex-end;
  gap: 26px;
  flex-wrap: wrap;
  justify-content: center;
}
.field {
  position: relative;
  width: clamp(180px, 17vw, 250px);
}
.field input {
  width: 100%;
  box-sizing: border-box;
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid #e6e6e4;
  padding: 10px 2px;
  font-size: 15.5px;
  font-family: inherit;
  color: #1a1a1a;
  letter-spacing: 0.02em;
  border-radius: 0;
}
.field input::placeholder {
  color: #cfcfcd;
  letter-spacing: 0.08em;
}
.line {
  position: absolute;
  left: 0;
  right: 0;
  bottom: 0;
  height: 1px;
  background: #1a1a1a;
  transform: scaleX(0);
  transition: transform 0.45s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.field:focus-within .line {
  transform: scaleX(1);
}
.cross {
  padding-bottom: 12px;
  color: #d7d7d5;
  font-size: 15px;
}
.go {
  width: 46px;
  height: 46px;
  margin-bottom: 1px;
  border: none;
  border-radius: 50%;
  background: #1a1a1a;
  color: #fdfdfd;
  display: grid;
  place-items: center;
  cursor: pointer;
  font-family: inherit;
  transition: transform 0.3s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.go:hover:not(:disabled) {
  transform: scale(1.06);
}
.go:disabled {
  opacity: 0.55;
  cursor: default;
}
.spin {
  width: 14px;
  height: 14px;
  border: 1.6px solid rgba(253, 253, 253, 0.35);
  border-top-color: #fdfdfd;
  border-radius: 50%;
  animation: rot 0.8s linear infinite;
}
@keyframes rot {
  to {
    transform: rotate(360deg);
  }
}
.retrieve {
  width: min(560px, 90vw);
  margin-top: 34px;
}
.failed-msg {
  margin-left: auto;
  margin-right: auto;
  font-size: 12px;
  letter-spacing: 0.05em;
  color: #b3261e;
  text-align: center;
}
.done {
  text-align: center;
}
.done-line {
  margin: 0 0 14px;
  font-size: 13.5px;
  letter-spacing: 0.06em;
  color: #1a1a1a;
}
.done-line b {
  font-weight: 600;
  font-variant-numeric: tabular-nums;
}
.failed {
  font-size: 11px;
  color: #b3261e;
}
.enter {
  margin: 0 auto;
}
.venn {
  max-width: min(520px, 90%);
  margin: 6px auto 0;
  display: block;
  background: #ffffff;
  border: 1px solid #f0f0ee;
}
.foot-stats {
  margin-top: 48px;
  font-size: 11px;
  letter-spacing: 0.1em;
  color: #bcbcbc;
}
.foot-stats b {
  font-weight: 500;
  color: #a3a3a3;
  font-variant-numeric: tabular-nums;
}
.hint-line {
  margin: 10px 0 0;
  font-size: 11px;
  letter-spacing: 0.06em;
  color: #c6c6c4;
}
</style>
