<script setup lang="ts">
/** 应用骨架:顶栏导航 / 入场幕布 / 呼吸细环 / 页脚(家族设计语言,样式见 styles/screen.css)
 * + hash 路由(未知 hash 走兜底页,不再空白)。 */
import { computed, onBeforeUnmount, onMounted, ref, watchEffect } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from './api'
import { setLang } from './i18n'
import { applyScale, persistScale, storedScale } from './utils/fontScale'
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

/* 下拉导航:分组名/条目名复用旧版 popover 词典(app.*),分组方式与 Streamlit 版一致。
 * 收藏 / ML 结果 / 数据管理三个从属页退出顶栏,入口移到各自宿主页(Docking/Ml/Analysis)。 */
interface NavItem { route: string; label: string }
interface NavGroup { key: string; label: string; items: NavItem[] }

const NAV_GROUPS: NavGroup[] = [
  {
    key: 'source', label: 'app.popover_source', items: [
      { route: 'drug-sources', label: 'app.btn_drug_mgmt' },
      { route: 'disease-sources', label: 'app.btn_disease_mgmt' },
      { route: 'ppi', label: 'app.btn_ppi_info' },
      { route: 'enrichment', label: 'app.btn_enrichment_info' }
    ]
  },
  {
    key: 'validation', label: 'app.popover_docking', items: [
      { route: 'analysis', label: 'app.btn_data_network' },
      { route: 'docking', label: 'app.btn_molecular_docking' },
      { route: 'ml', label: 'app.btn_expression_ml' }
    ]
  }
]
const NAV_DIRECT: NavItem[] = [
  { route: 'logs', label: 'web.nav.logs' },
  { route: 'settings', label: 'web.nav.settings' }
]

const lang = computed(() => (locale.value === 'zh' ? 'EN' : '中文'))
function toggleLang(): void {
  setLang(locale.value === 'zh' ? 'en' : 'zh')
}

/* <html lang> 跟随语言,便于屏幕阅读器与字体渲染 */
watchEffect(() => {
  document.documentElement.lang = locale.value === 'zh' ? 'zh-CN' : 'en'
})

const route = ref(location.hash.replace(/^#\/?/, '') || 'home')
const introVisible = ref(true)
const ready = ref(false)
let introTimer = 0

/* 下拉开合:点击开合;点击面板外 / Esc / 跳转后关闭 */
const openMenu = ref('')
function toggleMenu(key: string): void {
  openMenu.value = openMenu.value === key ? '' : key
}
function closeMenu(): void {
  openMenu.value = ''
}
function groupActive(g: NavGroup): boolean {
  return g.items.some((it) => it.route === route.value)
}

function onDocClick(e: MouseEvent): void {
  if (!openMenu.value) return
  const el = e.target instanceof Element ? e.target.closest('.dd') : null
  if (!el) closeMenu()
}
function onKeydown(e: KeyboardEvent): void {
  if (e.key === 'Escape') closeMenu()
}

function onHash(): void {
  route.value = location.hash.replace(/^#\/?/, '') || 'home'
  closeMenu()
}

onMounted(() => {
  window.addEventListener('hashchange', onHash)
  window.addEventListener('click', onDocClick)
  window.addEventListener('keydown', onKeydown)
  introTimer = window.setTimeout(() => {
    introVisible.value = false
    ready.value = true
  }, 1750)
  /* 界面缩放:localStorage 即时恢复,再以后端设置快照为准校准(另一浏览器保存过的值) */
  applyScale(storedScale())
  void api.settings()
    .then((s) => { persistScale(s.ui_font_scale); applyScale(s.ui_font_scale) })
    .catch(() => { /* 读取失败维持本地镜像值 */ })
})
onBeforeUnmount(() => {
  window.removeEventListener('hashchange', onHash)
  window.removeEventListener('click', onDocClick)
  window.removeEventListener('keydown', onKeydown)
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
        <a href="#/home" :class="{ active: route === 'home' }">{{ t('web.nav.home') }}</a>
        <div v-for="g in NAV_GROUPS" :key="g.key" class="dd" :class="{ on: groupActive(g) }">
          <button class="dd-btn" type="button" :aria-expanded="openMenu === g.key"
                  @click.stop="toggleMenu(g.key)">
            {{ t(g.label) }}<i class="dd-caret" aria-hidden="true"></i>
          </button>
          <Transition name="dd">
            <div v-if="openMenu === g.key" class="dd-panel">
              <a v-for="it in g.items" :key="it.route" class="dd-item"
                 :class="{ active: route === it.route }" :href="`#/${it.route}`"
                 @click="closeMenu">{{ t(it.label) }}</a>
            </div>
          </Transition>
        </div>
        <a v-for="d in NAV_DIRECT" :key="d.route" :href="`#/${d.route}`"
           :class="{ active: route === d.route }">{{ t(d.label) }}</a>
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
      <CollectionView v-else-if="route === 'ppi'" kind="ppi" :kicker="t('web.collection.ppi_kicker')"
                      :title="t('web.collection.ppi_title')" :sub="t('web.collection.ppi_sub')" />
      <CollectionView v-else-if="route === 'enrichment'" kind="enrichment"
                      :kicker="t('web.collection.enr_kicker')" :title="t('web.collection.enr_title')"
                      :sub="t('web.collection.enr_sub')" />
      <CollectionView v-else-if="route === 'stars'" kind="stars" :kicker="t('web.collection.stars_kicker')"
                      :title="t('web.collection.stars_title')" :sub="t('web.collection.stars_sub')" />
      <DataView v-else-if="route === 'data'" />
      <MlView v-else-if="route === 'ml'" />
      <MlResultsView v-else-if="route === 'ml-results'" />
      <SettingsView v-else-if="route === 'settings'" />
      <!-- 未知 hash 兜底:书签过期 / 手改地址时不至于白屏 -->
      <div v-else class="notfound">
        <p class="kicker s1">{{ t('web.notfound.kicker') }}</p>
        <h1 class="title s2">{{ t('web.notfound.title') }}</h1>
        <p class="sub s3">{{ t('web.notfound.sub', { hash: `#/${route}` }) }}</p>
        <a class="bar-btn nf-back s4" href="#/home">{{ t('web.notfound.back') }}</a>
      </div>
    </main>

    <footer class="foot">
      <span>{{ t('web.app.foot') }}</span>
    </footer>
  </div>
</template>

<style scoped>
.nf-back {
  margin-top: 34px;
}
</style>
