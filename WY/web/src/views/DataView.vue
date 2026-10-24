<script setup lang="ts">
/** 数据管理(阶段 2):背景互作组状态与一键重建(job)+ 分析产物目录(CollectionView)。 */
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { api, type InteractomeStatus } from '../api'
import { useJob } from '../composables/useJob'
import CollectionView from './CollectionView.vue'

const itx = ref<InteractomeStatus | null>(null)
const loaded = ref(false)
const notice = ref('')
let noticeTimer = 0

function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 3200)
}

const job = useJob((result) => {
  flash(`背景互作组构建完成:${result.nodes} 节点 / ${result.edges} 边` +
    ((result.cleaned as number) ? `,清理旧版本缓存 ${result.cleaned} 项` : ''))
  void load()
})

async function load(): Promise<void> {
  try {
    itx.value = await api.interactome()
    loaded.value = true
  } catch {
    flash('互作组状态读取失败')
  }
}

function build(): void {
  void job.run(api.buildInteractome(), { doneText: '构建完成' })
}

let pollTimer = 0
onMounted(() => {
  void load()
})
onBeforeUnmount(() => {
  window.clearTimeout(noticeTimer)
  if (!job.running.value) window.clearInterval(pollTimer)
})
</script>

<template>
  <p class="kicker s1">DATA MANAGEMENT</p>
  <h1 class="title s2">数据管理</h1>
  <p class="sub s3">背景互作组 · 分析产物目录</p>

  <div class="panel s4">
    <div class="controls">
      <span class="ctl-label">背景互作组</span>
      <span v-if="loaded && itx?.built" class="count">
        {{ itx.nodes }} 节点 · {{ itx.edges }} 边 · STRING {{ itx.version }} · 置信度
        {{ itx.confidence }} · 构建于 {{ itx.built_at }}(已 {{ itx.age ?? '-' }} 天)
      </span>
      <span v-else-if="loaded" class="count" style="color:#b3261e">
        未构建 —— 邻近度与 RWR 需要它,请先构建
      </span>
      <button class="go" title="重建" :disabled="job.running.value" @click="build">
        <svg viewBox="0 0 20 20" width="15" height="15">
          <path d="M16.5 10a6.5 6.5 0 1 1-1.9-4.6M16.5 2.5v3.4h-3.4" fill="none"
                stroke="currentColor" stroke-width="1.6" stroke-linecap="round"
                stroke-linejoin="round" />
        </svg>
      </button>
    </div>

    <p v-if="loaded && itx?.built && itx?.stale" class="notice">
      背景互作组已 {{ itx.age }} 天,数据可能陈旧 —— 建议重建。
    </p>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ job.text.value || '构建中…(首次需下载数 GB,可能数十分钟)' }}</p>
    </div>
    <p v-if="job.error.value" class="notice">构建失败:{{ job.error.value }}</p>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>

  <CollectionView kind="downloads" kicker="ANALYSIS ARTIFACTS" title="分析产物目录"
                  sub="每次分析归档一个目录 · 分子对接成果与收藏由各自页面管理" />
</template>
