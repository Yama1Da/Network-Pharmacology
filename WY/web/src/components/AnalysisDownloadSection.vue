<script setup lang="ts">
/** 分析页下载区:运行目录全量文件清单(全量枚举保底 + 友好名叠加)+ 整目录 zip。
 * 纯展示 + 单一动作;清单数据由父级从 ppi.files 下发(一致性面板 @saved 会追加
 * 文件进父级清单,经 props 传导到这里)。样式 .fold/.fold-title 已上提全局。 */
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, artifactsUrl } from '../api'
import { downloadBlob } from '../utils/download'
import { dlDesc } from '../utils/dlFiles'

const props = defineProps<{
  files: { name: string; rel: string }[]
  directoryRel: string
  flash: (msg: string) => void
}>()
const { t } = useI18n()

const downloadItems = computed(() =>
  props.files.map((f) => ({
    name: f.name, rel: f.rel, desc: dlDesc(f.name, t),
    known: dlDesc(f.name, t) !== f.name,
  })))

const dirName = computed(() =>
  (props.directoryRel ?? '').split('/').pop() ?? '')

async function zipWholeDirectory(): Promise<void> {
  if (!dirName.value) return
  try {
    downloadBlob(await api.zip('downloads', [dirName.value]), `${dirName.value}.zip`)
  } catch {
    props.flash(t('web.analysis.dl_zip_fail'))
  }
}
</script>

<template>
  <section class="sec">
    <details class="fold">
      <summary class="sec-title fold-title">{{ t('web.analysis.dl_title', { n: downloadItems.length }) }}</summary>
      <button class="bar-btn" @click="zipWholeDirectory">{{ t('web.analysis.dl_whole', { name: `${dirName}.zip` }) }}</button>
      <ul class="dl-list">
        <li v-for="d in downloadItems" :key="d.rel" class="dl-row">
          <div class="dl-main">
            <span class="dl-desc">{{ d.desc }}</span>
            <span class="dl-name" :class="{ unknown: !d.known }">{{ d.name }}</span>
          </div>
          <a class="mini" :href="artifactsUrl(d.rel, true)">{{ t('web.common.download') }}</a>
        </li>
      </ul>
      <p class="hint-line">{{ t('web.analysis.dl_hint', { dir: directoryRel }) }}</p>
    </details>
  </section>
</template>
