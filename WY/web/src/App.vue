<script setup lang="ts">
/** 应用骨架:顶栏导航 / 入场幕布 / 呼吸细环 / 页脚(家族设计语言,样式见 styles/screen.css)
 * + hash 路由。过渡期 Streamlit(8501)仍是主入口,这里承载逐页迁入的 TS 页面。 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { setLang } from './i18n'
import CollectionView from './views/CollectionView.vue'
import AnalysisView from './views/AnalysisView.vue'
import DataView from './views/DataView.vue'
import DockingView from './views/DockingView.vue'
import HomeView from './views/HomeView.vue'
import MlView from './views/MlView.vue'
import MlResultsView from './views/MlResultsView.vue'
import LogsView from './views/LogsView.vue'
import SettingsView from './views/SettingsView.vue'
import SourcesView from './views/SourcesView.vue'

const { t, locale } = useI18n()

const NAV = [
  { hash: '#/home', key: 'web.nav.home' },
  { hash: '#/analysis', key: 'web.nav.analysis' },
  { hash: '#/docking', key: 'web.nav.docking' },
  { hash: '#/data', key: 'web.nav.data' },
  { hash: '#/logs', key: 'web.nav.logs' },
  { hash: '#/drug-sources', key: 'web.nav.drug_sources' },
  { hash: '#/disease-sources', key: 'web.nav.disease_sources' },
  { hash: '#/ppi', key: 'web.nav.ppi' },
  { hash: '#/enrichment', key: 'web.nav.enrichment' },
  { hash: '#/stars', key: 'web.nav.stars' },
  { hash: '#/ml', key: 'web.nav.ml' },
  { hash: '#/ml-results', key: 'web.nav.ml_results' },
  { hash: '#/settings', key: 'web.nav.settings' },
]

const lang = computed(() => (locale.value === 'zh' ? 'EN' : '中文'))
function toggleLang(): void {
  setLang(locale.value === 'zh' ? 'en' : 'zh')
}

const route = ref(location.hash.replace(/^#\/?/, '') || 'home')
const introVisible = ref(true)
const ready = ref(false)
let introTimer = 0

function onHash(): void {
  route.value = location.hash.replace(/^#\/?/, '') || 'home'
}

onMounted(() => {
  window.addEventListener('hashchange', onHash)
  introTimer = window.setTimeout(() => {
    introVisible.value = false
    ready.value = true
  }, 1750)
})
onBeforeUnmount(() => {
  window.removeEventListener('hashchange', onHash)
  window.clearTimeout(introTimer)
})
</script>

<template>
  <div class="screen">
    <div class="whisper w2" aria-hidden="true"></div>
    <div class="whisper w1" aria-hidden="true"><i class="orbiter"></i></div>

    <Transition name="intro-fade">
      <div v-if="introVisible" class="intro" aria-hidden="true">
        <svg class="intro-ring" viewBox="0 0 64 64" width="56" height="56">
          <circle class="intro-ring-path" cx="32" cy="32" r="30" />
        </svg>
        <span class="intro-word">NPAP</span>
      </div>
    </Transition>

    <header class="top">
      <a class="brand" href="#/home">NPAP<i class="brand-dot"></i></a>
      <nav class="links">
        <a v-for="n in NAV" :key="n.hash" :href="n.hash"
           :class="{ active: route === n.hash.slice(2) }">{{ t(n.key) }}</a>
        <a @click.prevent="toggleLang">{{ lang }}</a>
      </nav>
    </header>

    <main class="center" :class="{ ready }">
      <HomeView v-if="route === 'home'" />
      <AnalysisView v-else-if="route === 'analysis'" />
      <DockingView v-else-if="route === 'docking'" />
      <LogsView v-else-if="route === 'logs'" />
      <SourcesView v-else-if="route === 'drug-sources'" kind="drug" />
      <SourcesView v-else-if="route === 'disease-sources'" kind="disease" />
      <CollectionView v-else-if="route === 'ppi'" kind="ppi" kicker="PPI NETWORK FILES"
                      title="PPI 文件" sub="STRING 互作缓存 · 可打包可删除" />
      <CollectionView v-else-if="route === 'enrichment'" kind="enrichment"
                      kicker="ENRICHMENT CACHE" title="富集文件"
                      sub="g:Profiler · Enrichr · GSEA 缓存" />
      <CollectionView v-else-if="route === 'stars'" kind="stars" kicker="STARRED RESULTS"
                      title="收藏的对接结果" sub="分子对接收藏 · 目录整体打包" />
      <DataView v-else-if="route === 'data'" />
      <MlView v-else-if="route === 'ml'" />
      <MlResultsView v-else-if="route === 'ml-results'" />
      <SettingsView v-else-if="route === 'settings'" />
    </main>

    <footer class="foot">
      <span>© 2026 NPAP · 仅供科研使用</span>
    </footer>
  </div>
</template>
