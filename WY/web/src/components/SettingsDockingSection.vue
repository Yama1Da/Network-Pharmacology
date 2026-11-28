<script setup lang="ts">
/** 设置页 ③ 分子对接:口袋阈值/Vina 参数(走 /api/settings/ui)+ Vina·Open Babel
 * 路径(走 save_api_env)双通道保存。快照由壳层下发(watch immediate 回填本栏
 * ref);保存成功 emit('saved', 新快照) 由壳层统一回填各栏。 */
import { ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot } from '../api'
import Stepper from './Stepper.vue'
import { useSaveNote } from '../composables/useSaveNote'

const props = defineProps<{ snapshot: SettingsSnapshot | null }>()
const emit = defineEmits<{ saved: [SettingsSnapshot] }>()
const { t } = useI18n()
const { saveNote, noteFlash, flashSaved } = useSaveNote()

const pocketCutoff = ref(4.0)
const exhaustiveness = ref(8)
const dockingCpu = ref(4)
const vinaDir = ref('')
const obabelDir = ref('')

watch(() => props.snapshot, (s) => {
  if (!s) return
  pocketCutoff.value = s.pocket_cutoff
  exhaustiveness.value = s.exhaustiveness
  dockingCpu.value = s.docking_cpu
  vinaDir.value = s.vina_dir
  obabelDir.value = s.obabel_dir
}, { immediate: true })

async function saveDocking(): Promise<void> {
  try {
    const ui = await api.saveUiSettings({
      POCKET_CUTOFF: pocketCutoff.value,
      EXHAUSTIVENESS: exhaustiveness.value,
      DOCKING_CPU: dockingCpu.value
    })
    emit('saved', ui.snapshot)
    const applied = (await api.saveApiEnv({
      VINA_DIR: vinaDir.value, OBABEL_DIR: obabelDir.value
    })).applied
    flashSaved(0)
    if (applied.length) noteFlash(true, `${t('web.settings.saved')} · ${applied.join(', ')}`)
  } catch (e) {
    noteFlash(false, e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}
</script>

<template>
  <div class="opt-grid">
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.pocket_label') }}</span>
      <Stepper v-model="pocketCutoff" :min="3" :max="6" :step="0.5" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.exhaust_label') }}</span>
      <Stepper v-model="exhaustiveness" :min="1" :max="32" :step="1" :digits="0" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.cpu_label') }}</span>
      <Stepper v-model="dockingCpu" :min="0" :max="64" :step="1" :digits="0" />
    </label>
  </div>
  <p class="hint-line">{{ t('web.settings.exhaust_hint') }}</p>
  <p class="hint-line">{{ t('web.settings.cpu_hint') }}</p>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.sec_paths') }}</h3>
    <p class="hint-line">{{ t('web.settings.paths_hint') }}</p>
    <div class="set-row">
      <label class="set-label">{{ t('web.settings.vina_label') }}</label>
      <div class="field wide"><input v-model="vinaDir" spellcheck="false" /><i class="line"></i></div>
    </div>
    <div class="set-row">
      <label class="set-label">{{ t('web.settings.obabel_label') }}</label>
      <div class="field wide"><input v-model="obabelDir" spellcheck="false" /><i class="line"></i></div>
    </div>
  </section>
  <button class="bar-btn" @click="saveDocking">{{ t('web.settings.save') }}</button>
  <Transition name="note">
    <span v-if="saveNote" class="save-note" :class="{ fail: !saveNote.ok }">{{ saveNote.text }}</span>
  </Transition>
</template>
