<script setup lang="ts">
/** 设置页 ④ 通用设置:API 密钥/CLI 产物目录(save_api_env 双通道)+ 字号/主题
 * 显示偏好(即时生效 + 本地镜像,不触发缓存清理);含数据迁移(SettingsTransfer)
 * 挂载与导入界面设置的回填。快照由壳层下发(watch immediate 回填本栏 ref);
 * 保存成功 emit('saved')(无载荷)由壳层重新拉取快照(save_api_env 不返回快照)。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot } from '../api'
import SettingsTransfer from './SettingsTransfer.vue'
import { useSaveNote } from '../composables/useSaveNote'
import { applyScale, persistScale, snapScale } from '../utils/fontScale'
import { normalizeLang, setLang } from '../i18n'
import { applyTheme, normalizeTheme, persistTheme, type Theme } from '../utils/theme'

const props = defineProps<{
  snapshot: SettingsSnapshot | null
  /** 顶部 notice(SettingsTransfer 的导出/导入流程走顶部提示) */
  flash: (msg: string) => void
  /** 迁移落盘后的全量刷新 */
  reload: () => Promise<void>
}>()
const emit = defineEmits<{ saved: [] }>()
const { t } = useI18n()
const { saveNote, noteFlash } = useSaveNote()

const keys = ref({ disgenet: false, ncbi: false })
const dgKey = ref('')
const ncbiKey = ref('')
const cliResultsDir = ref('')

/* 字体大小(界面等比缩放):滑动条 80%–150%、步长 10%,savedScale 用于保存按钮的未修改判定 */
const fontScale = ref(1)
const savedScale = ref(1)
const fontScaleDirty = computed(() => fontScale.value !== savedScale.value)

/* 颜色变换:浅色/深色(样品「墨黑」),savedTheme 用于保存按钮的未修改判定 */
const themeSel = ref<Theme>('light')
const savedTheme = ref<Theme>('light')
const themeDirty = computed(() => themeSel.value !== savedTheme.value)

watch(() => props.snapshot, (s) => {
  if (!s) return
  keys.value = s.keys
  cliResultsDir.value = s.cli_results_dir
  fontScale.value = snapScale(s.ui_font_scale)
  savedScale.value = fontScale.value
  themeSel.value = normalizeTheme(s.ui_theme)
  savedTheme.value = themeSel.value
}, { immediate: true })

async function saveGeneral(): Promise<void> {
  try {
    if (fontScaleDirty.value || themeDirty.value) {
      // 缩放/主题即时生效 + 本地镜像,再随快照确认后端值(纯显示偏好,不触发缓存清理)
      const ui = await api.saveUiSettings({
        ...(fontScaleDirty.value ? { UI_FONT_SCALE: fontScale.value } : {}),
        ...(themeDirty.value ? { UI_THEME: themeSel.value } : {})
      })
      persistScale(fontScale.value)
      applyScale(fontScale.value)
      persistTheme(themeSel.value)
      applyTheme(themeSel.value)
    }
    const changes: Record<string, string> = { CLI_RESULTS_DIR: cliResultsDir.value }
    if (dgKey.value.trim()) changes.DisGeNET_API_KEY = dgKey.value
    if (ncbiKey.value.trim()) changes.NCBI_API_KEY = ncbiKey.value
    const applied = (await api.saveApiEnv(changes)).applied
    dgKey.value = ''
    ncbiKey.value = ''
    noteFlash(true, applied.length ? t('web.settings.keys_saved', { names: applied.join(', ') })
                                   : t('web.settings.saved'))
    emit('saved')
  } catch {
    noteFlash(false, t('web.settings.keys_fail'))
  }
}

/* 导入的界面设置(主题/字号/语言)回填:后端已生效,前端同步本地镜像与当前显示 */
function applyImportedUi(snap: SettingsSnapshot): void {
  themeSel.value = normalizeTheme(snap.ui_theme)
  persistTheme(themeSel.value)
  applyTheme(themeSel.value)
  fontScale.value = snapScale(snap.ui_font_scale)
  savedScale.value = fontScale.value
  persistScale(fontScale.value)
  applyScale(fontScale.value)
  setLang(normalizeLang(snap.ui_lang))
}
</script>

<template>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.sec_keys') }}</h3>
    <p class="hint-line">{{ t('web.settings.keys_hint') }}</p>
    <div class="set-row">
      <label class="set-label">DisGeNET API Key</label>
      <div class="field"><input v-model="dgKey" type="password" placeholder="********"
                                spellcheck="false" /><i class="line"></i></div>
      <span class="state" :class="{ on: keys.disgenet }">{{ keys.disgenet ? t('web.settings.configured') : t('web.settings.not_configured') }}</span>
    </div>
    <div class="set-row">
      <label class="set-label">NCBI API Key</label>
      <div class="field"><input v-model="ncbiKey" type="password" placeholder="********"
                                spellcheck="false" /><i class="line"></i></div>
      <span class="state" :class="{ on: keys.ncbi }">{{ keys.ncbi ? t('web.settings.configured') : t('web.settings.not_configured') }}</span>
    </div>
  </section>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.cli_title') }}</h3>
    <div class="set-row">
      <label class="set-label">{{ t('web.settings.cli_label') }}</label>
      <div class="field wide"><input v-model="cliResultsDir" spellcheck="false" /><i class="line"></i></div>
    </div>
  </section>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.font_title') }}</h3>
    <p class="hint-line">{{ t('web.settings.font_hint') }}</p>
    <div class="range-row">
      <input v-model.number="fontScale" type="range" min="0.8" max="1.5" step="0.1"
             :aria-label="t('web.settings.font_title')" />
      <span class="range-val">{{ Math.round(fontScale * 100) }}%</span>
    </div>
  </section>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.theme_title') }}</h3>
    <p class="hint-line">{{ t('web.settings.theme_hint') }}</p>
    <div class="ver-row">
      <button class="lv" :class="{ on: themeSel === 'light' }"
              @click="themeSel = 'light'">{{ t('web.settings.theme_light') }}</button>
      <span class="lv-sep">/</span>
      <button class="lv" :class="{ on: themeSel === 'dark' }"
              @click="themeSel = 'dark'">{{ t('web.settings.theme_dark') }}</button>
    </div>
  </section>
  <button class="bar-btn"
          :disabled="!dgKey.trim() && !ncbiKey.trim() && cliResultsDir === '' && !fontScaleDirty && !themeDirty"
          @click="saveGeneral">{{ t('web.settings.save') }}</button>
  <Transition name="note">
    <span v-if="saveNote" class="save-note" :class="{ fail: !saveNote.ok }">{{ saveNote.text }}</span>
  </Transition>

  <SettingsTransfer :flash="flash" :reload="reload" @ui-imported="applyImportedUi" />
</template>
