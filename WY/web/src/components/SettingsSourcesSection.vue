<script setup lang="ts">
/** 设置页 ① 数据来源:抓取规模(OT/DisGeNET 页数)与 STRING 版本切换。
 * 快照由壳层下发(watch immediate 回填本栏 ref);保存成功 emit('saved', 新快照)
 * 由壳层统一回填各栏——对应拆分前 applySnapshot(r.snapshot) 的全量回填语义。 */
import { ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot } from '../api'
import Stepper from './Stepper.vue'
import { useSaveNote } from '../composables/useSaveNote'

const props = defineProps<{ snapshot: SettingsSnapshot | null }>()
const emit = defineEmits<{ saved: [SettingsSnapshot] }>()
const { t } = useI18n()
const { saveNote, noteFlash, flashSaved } = useSaveNote()

const otTopN = ref(500)
const disgenetPages = ref(5)
const stringVersion = ref('')
const knownVersions = ref<string[]>([])
const customMode = ref(false)
const customVersion = ref('')

watch(() => props.snapshot, (s) => {
  if (!s) return
  otTopN.value = s.ot_top_n
  disgenetPages.value = s.disgenet_pages
  stringVersion.value = s.string_version
  knownVersions.value = s.known_versions
}, { immediate: true })

async function saveSources(): Promise<void> {
  try {
    const r = await api.saveUiSettings({
      OT_TOP_N: otTopN.value, DISGENET_PAGES: disgenetPages.value
    })
    emit('saved', r.snapshot)
    flashSaved(r.purged_merged)
  } catch (e) {
    noteFlash(false, e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}

async function applyStringVersion(): Promise<void> {
  const target = customMode.value ? customVersion.value.trim() : stringVersion.value
  if (!target) return
  try {
    const r = await api.setStringVersion(target)
    noteFlash(true, r.changed ? t('web.settings.switched', { v: r.version }) : t('web.settings.no_change'))
    stringVersion.value = r.version
  } catch (e) {
    noteFlash(false, e instanceof Error ? e.message : t('web.settings.switch_fail'))
  }
}

function chooseVersion(v: string): void {
  stringVersion.value = v
}
</script>

<template>
  <div class="opt-grid">
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.ot_top_n') }}</span>
      <Stepper v-model="otTopN" :min="200" :max="1000" :step="50" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.disgenet_pages') }}</span>
      <Stepper v-model="disgenetPages" :min="1" :max="10" :step="1" />
    </label>
  </div>

  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.sec_string', { v: stringVersion }) }}</h3>
    <p class="hint-line">{{ t('web.settings.string_hint') }}</p>
    <div class="ver-row">
      <button v-for="v in knownVersions" :key="v" class="lv" :class="{ on: stringVersion === v }"
              @click="chooseVersion(v)">{{ v }}</button>
      <button class="lv" :class="{ on: customMode }" @click="customMode = !customMode">{{ t('web.settings.custom_btn') }}</button>
    </div>
    <div v-if="customMode" class="set-row">
      <label class="set-label">{{ t('web.settings.version_label') }}</label>
      <div class="field"><input v-model="customVersion" placeholder="v13.0" spellcheck="false" /><i class="line"></i></div>
      <button class="bar-btn" :disabled="!customVersion.trim()" @click="applyStringVersion">{{ t('web.settings.switch_custom') }}</button>
    </div>
    <div v-else-if="stringVersion && !knownVersions.includes(stringVersion)" class="ver-row">
      <span class="hint-line">{{ t('web.settings.unknown_version', { v: stringVersion }) }}</span>
    </div>
  </section>

  <button class="bar-btn" @click="saveSources">{{ t('web.settings.save') }}</button>
  <Transition name="note">
    <span v-if="saveNote" class="save-note" :class="{ fail: !saveNote.ok }">{{ saveNote.text }}</span>
  </Transition>
</template>
