<script setup lang="ts">
/** 分析页网络图:视图切换(interactive iframe / static / STRING / Cytoscape)
 * + 核心靶点图。纯展示,views/artifacts 由父级从 PpiResult 下发,产物经
 * /api/artifacts 伺服。默认视图:有 cytoscape 版选 cytoscape,否则 interactive。 */
import { computed, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { artifactsUrl } from '../api'

const props = defineProps<{
  views: string[]
  artifacts: Record<string, string>
}>()
const { t } = useI18n()

const view = ref('')
const effectiveView = computed(() => {
  const v = view.value
  if (v && props.views.includes(v)) return v
  return props.views.includes('cytoscape') ? 'cytoscape' : 'interactive'
})
const viewLabels = computed<Record<string, string>>(() => ({
  interactive: t('web.analysis.view_interactive'),
  static: t('web.analysis.view_static'),
  string: t('web.analysis.view_string'),
  cytoscape: t('web.analysis.view_cytoscape')
}))
function art(key: string): string | null {
  return props.artifacts[key] ?? null
}
</script>

<template>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.analysis.sec_network') }}</h3>
    <div class="lv-group">
      <button v-for="v in views" :key="v" class="lv"
              :class="{ on: effectiveView === v }" @click="view = v">{{ viewLabels[v] }}</button>
    </div>
    <p v-if="effectiveView === 'string'" class="hint-line">{{ t('web.analysis.view_string_hint') }}</p>
    <p v-if="effectiveView === 'cytoscape'" class="hint-line">{{ t('web.analysis.view_cyto_hint') }}</p>
    <iframe v-if="effectiveView === 'interactive'"
            class="fig-frame" :src="artifactsUrl(art('html') ?? '')" />
    <img v-else-if="effectiveView === 'static'"
         class="fig-img" :src="artifactsUrl(art('png') ?? '')" alt="static network" />
    <img v-else-if="effectiveView === 'string'"
         class="fig-img" :src="artifactsUrl(art('string_png') ?? '')" alt="string style" />
    <img v-else class="fig-img" :src="artifactsUrl(art('cyto_png') ?? '')" alt="cytoscape style" />
    <p class="hint-line">{{ t('web.analysis.views_saved_hint') }}</p>
  </section>

  <section v-if="art('core_png')" class="sec">
    <h3 class="sec-title">{{ t('web.analysis.sec_core') }}</h3>
    <img class="fig-img" :src="artifactsUrl(art('core_png') ?? '')" alt="core targets" />
  </section>
</template>
