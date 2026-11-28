<script setup lang="ts">
/** 设置页 ② 分析算法:ChEMBL 黑名单、降噪/富集阈值、统计稳健性参数(缓存键均含
 * 这些值,改动后旧缓存自动重算)。保存走 /api/settings/ui(后端清 merged 整合缓存
 * 并同步互作组置信度)。快照由壳层下发(watch immediate 回填本栏 ref);保存成功
 * emit('saved', 新快照) 由壳层统一回填各栏。 */
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot } from '../api'
import Stepper from './Stepper.vue'
import { useSaveNote } from '../composables/useSaveNote'

const props = defineProps<{ snapshot: SettingsSnapshot | null }>()
const emit = defineEmits<{ saved: [SettingsSnapshot] }>()
const { t } = useI18n()
const { saveNote, noteFlash, flashSaved } = useSaveNote()

const blacklist = ref<string[]>([])
const blacklistOptions = ref<string[]>([])
const newGene = ref('')
const pchembl = ref(0.5)
const dgidbMin = ref(0.05)
const corroboration = ref(false)
const otMin = ref(0.1)
const enrichP = ref(0.055)
const louvain = ref(1.0)
const ppiConf = ref(0.4)

/* 统计稳健性/研究语义参数 */
const rwrAlpha = ref(0.7)
const gseaPerms = ref(1000)
const proxPerms = ref(200)
const rwrNullPerms = ref(200)
const synergyTopN = ref(100)
const synergyPerms = ref(200)
const gseaMinSize = ref(15)
const gseaMaxSize = ref(500)
const gseaFdr = ref(0.25)
const enrichSources = ref<string[]>([])
const ENRICHMENT_SOURCE_VALUES = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG']
/* 交叉校验:规模窗不允许倒挂、富集类别不允许空集(后端同款校验,前端先行拦截) */
const gseaSizeBad = computed(() => gseaMinSize.value >= gseaMaxSize.value)
const sourcesEmpty = computed(() => enrichSources.value.length === 0)

watch(() => props.snapshot, (s) => {
  if (!s) return
  blacklist.value = s.chembl_blacklist
  blacklistOptions.value = s.chembl_blacklist_options
  pchembl.value = s.chembl_pchembl_threshold
  dgidbMin.value = s.dgidb_min_score
  corroboration.value = s.chembl_corroboration_required
  otMin.value = s.ot_min_score
  enrichP.value = s.enrichment_pvalue
  louvain.value = s.louvain_resolution
  ppiConf.value = s.ppi_confidence
  rwrAlpha.value = s.rwr_alpha
  gseaPerms.value = s.gsea_permutation_num
  proxPerms.value = s.proximity_n_perm
  rwrNullPerms.value = s.rwr_null_perms
  synergyTopN.value = s.synergy_top_n
  synergyPerms.value = s.synergy_n_perm
  gseaMinSize.value = s.gsea_set_min_size
  gseaMaxSize.value = s.gsea_set_max_size
  gseaFdr.value = s.gsea_fdr_screen
  enrichSources.value = [...s.enrichment_sources]
}, { immediate: true })

const candidateOptions = computed(() =>
  blacklistOptions.value.filter((o) => !blacklist.value.includes(o)))

function addGene(gene?: string): void {
  const g = (gene ?? newGene.value).trim().toUpperCase()
  if (g && !blacklist.value.includes(g)) blacklist.value = [...blacklist.value, g]
  newGene.value = ''
}

function removeGene(gene: string): void {
  blacklist.value = blacklist.value.filter((g) => g !== gene)
}

async function saveAlgo(): Promise<void> {
  if (sourcesEmpty.value) {
    noteFlash(false, t('web.settings.sources_empty'))
    return
  }
  if (gseaSizeBad.value) {
    noteFlash(false, t('web.settings.gsea_size_bad'))
    return
  }
  try {
    const r = await api.saveUiSettings({
      CHEMBL_TARGET_BLACKLIST: blacklist.value,
      CHEMBL_PCHEMBL_THRESHOLD: pchembl.value,
      DGIDB_MIN_SCORE: dgidbMin.value,
      CHEMBL_CORROBORATION_REQUIRED: corroboration.value,
      OT_MIN_SCORE: otMin.value,
      ENRICHMENT_PVALUE: enrichP.value,
      LOUVAIN_RESOLUTION: louvain.value,
      PPI_CONFIDENCE: ppiConf.value,
      RWR_ALPHA: rwrAlpha.value,
      GSEA_PERMUTATION_NUM: gseaPerms.value,
      PROXIMITY_N_PERM: proxPerms.value,
      RWR_NULL_PERMS: rwrNullPerms.value,
      SYNERGY_TOP_N: synergyTopN.value,
      SYNERGY_N_PERM: synergyPerms.value,
      GSEA_SET_MIN_SIZE: gseaMinSize.value,
      GSEA_SET_MAX_SIZE: gseaMaxSize.value,
      GSEA_FDR_SCREEN: gseaFdr.value,
      ENRICHMENT_SOURCES: enrichSources.value
    })
    emit('saved', r.snapshot)
    flashSaved(r.purged_merged)
  } catch (e) {
    noteFlash(false, e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}
</script>

<template>
  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.blacklist_label') }}</h3>
    <div v-if="blacklist.length" class="chips">
      <span v-for="g in blacklist" :key="g" class="chip">
        {{ g }}<button type="button" class="chip-x" :aria-label="`${t('web.settings.blacklist_remove')} ${g}`"
                       @click="removeGene(g)">×</button>
      </span>
    </div>
    <p v-else class="hint-line">{{ t('web.settings.blacklist_empty') }}</p>
    <div class="add-row">
      <div class="field">
        <input v-model="newGene" spellcheck="false" :placeholder="t('web.settings.blacklist_add_ph')"
               @keydown.enter.prevent="addGene()" />
        <i class="line"></i>
      </div>
      <button class="bar-btn" :disabled="!newGene.trim()" @click="addGene()">{{ t('web.settings.blacklist_add') }}</button>
    </div>
    <template v-if="candidateOptions.length">
      <p class="hint-line">{{ t('web.settings.blacklist_candidates') }}</p>
      <div class="chips">
        <button v-for="o in candidateOptions" :key="o" type="button" class="chip cand"
                @click="addGene(o)">{{ o }} +</button>
      </div>
    </template>
  </section>

  <div class="opt-grid">
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.pchembl_label') }}</span>
      <Stepper v-model="pchembl" :min="0" :max="1" :step="0.05" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.dgidb_label') }}</span>
      <Stepper v-model="dgidbMin" :min="0" :max="1" :step="0.01" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.ot_min_label') }}</span>
      <Stepper v-model="otMin" :min="0" :max="1" :step="0.01" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.enrich_p_label') }}</span>
      <Stepper v-model="enrichP" :min="0.01" :max="0.2" :step="0.005" :digits="3" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.louvain_label') }}</span>
      <Stepper v-model="louvain" :min="0.5" :max="2" :step="0.1" />
    </label>
    <label class="opt">
      <span class="opt-label">{{ t('web.settings.ppi_conf_label') }}</span>
      <Stepper v-model="ppiConf" :min="0.15" :max="0.9" :step="0.05" />
    </label>
  </div>
  <div class="corro-row">
    <label class="chk">
      <input v-model="corroboration" type="checkbox" /> {{ t('web.settings.corroboration_label') }}
    </label>
  </div>

  <section class="sec">
    <h3 class="sec-title">{{ t('web.settings.sec_stats') }}</h3>
    <p class="hint-line">{{ t('web.settings.stats_hint') }}</p>
    <div class="opt-grid">
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.rwr_alpha_label') }}</span>
        <Stepper v-model="rwrAlpha" :min="0.1" :max="0.9" :step="0.05" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.gsea_perms_label') }}</span>
        <Stepper v-model="gseaPerms" :min="100" :max="5000" :step="100" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.prox_perms_label') }}</span>
        <Stepper v-model="proxPerms" :min="50" :max="2000" :step="50" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.rwr_null_label') }}</span>
        <Stepper v-model="rwrNullPerms" :min="50" :max="2000" :step="50" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.synergy_top_n_label') }}</span>
        <Stepper v-model="synergyTopN" :min="5" :max="1000" :step="10" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.synergy_perms_label') }}</span>
        <Stepper v-model="synergyPerms" :min="10" :max="2000" :step="50" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.gsea_min_label') }}</span>
        <Stepper v-model="gseaMinSize" :min="5" :max="100" :step="5" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.gsea_max_label') }}</span>
        <Stepper v-model="gseaMaxSize" :min="100" :max="1000" :step="100" :digits="0" />
      </label>
      <label class="opt">
        <span class="opt-label">{{ t('web.settings.gsea_fdr_label') }}</span>
        <Stepper v-model="gseaFdr" :min="0.05" :max="0.5" :step="0.05" />
      </label>
    </div>
    <p class="hint-line">{{ t('web.settings.gsea_size_hint') }}</p>

    <p class="hint-line">{{ t('web.settings.sources_label') }}</p>
    <div class="tf-scopes">
      <label v-for="s in ENRICHMENT_SOURCE_VALUES" :key="s" class="chk">
        <input v-model="enrichSources" type="checkbox" :value="s" /> {{ s }}
      </label>
    </div>
  </section>

  <button class="bar-btn" :disabled="sourcesEmpty || gseaSizeBad"
          @click="saveAlgo">{{ t('web.settings.save') }}</button>
  <Transition name="note">
    <span v-if="saveNote" class="save-note" :class="{ fail: !saveNote.ok }">{{ saveNote.text }}</span>
  </Transition>
</template>
