<script setup lang="ts">
/** 协同预测页(二期):药物对表单 + 页级参数(N/置换轮数,Settings 作初始默认,页内覆盖)
 * + jobs 进度(stage 码映射文案,fraction 由核心按零分布轮数推进)+ 契约 synergy.v1.1 结果渲染。
 * 结果页联动:归档落 synergy-results/,「协同分析结果管理」提供历史清单/复现回填。
 * 措辞红线(§1.3):geometry 三枚举只述几何事实,禁译"协同/拮抗"。
 * 深色适配:仅家族类 + 布局 scoped 样式,色值全走 var(--c-*),无硬编码。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SynergyContract, type SynergyJobResult } from '../api'
import SynergyResultPanel from '../components/SynergyResultPanel.vue'
import { useJob } from '../composables/useJob'
import { useNotice } from '../composables/useNotice'
import { downloadBlob } from '../utils/download'
import { takePrefill, takeSynergyReplay } from '../store'

const { t } = useI18n()
const { notice, flash } = useNotice(3200)

const drugA = ref('')
const drugB = ref('')
const disease = ref('')
const drugOptions = ref<string[]>([])
const diseaseOptions = ref<string[]>([])
const N_OPTIONS: (number | null)[] = [30, 50, 100, null]      // null = 全量
const PERM_OPTIONS = [50, 100, 200, 400]
const topN = ref<number | null>(100)
const nPerm = ref<number | null>(200)
const nPermDefault = ref<number>(200)                          // Settings 初始默认,展示兜底用

const result = ref<SynergyContract | null>(null)
const archiveDir = ref('')

const job = useJob((r) => {
  const payload = r as unknown as SynergyJobResult
  result.value = payload
  archiveDir.value = payload.directory_rel ?? ''
})

onMounted(async () => {
  // 复现回填优先(结果管理页「复现此运行」),其次主页检索预填
  const replay = takeSynergyReplay()
  if (replay) {
    drugA.value = replay.drug_a
    drugB.value = replay.drug_b
    disease.value = replay.disease
    topN.value = replay.top_n
    nPerm.value = replay.n_perm ?? 200
  } else {
    const pf = takePrefill()
    if (pf) {
      drugA.value = pf.drug
      disease.value = pf.disease
    }
  }
  try {
    const lib = await api.synergyLibrary()
    drugOptions.value = lib.drugs
    diseaseOptions.value = lib.diseases
  } catch {
    /* 本地库清单不可达时输入框照常手输 */
  }
  try {
    const s = await api.settings()
    if (!replay) {
      topN.value = s.synergy_top_n ?? 100
      nPerm.value = s.synergy_n_perm ?? 200
    }
    nPermDefault.value = s.synergy_n_perm ?? 200
  } catch {
    /* 设置不可达时保持当前值 */
  }
  // 在途任务恢复(跨页离开任务照跑,重进页面找回进度与结果——同 PK 镜像语义)
  try {
    const st = await api.synergyJob()
    if (st.running && st.job_id) {
      void job.resume(st.job_id)
    }
  } catch {
    /* 无在途任务或状态不可达 */
  }
})

function run(): void {
  if (!drugA.value.trim() || !drugB.value.trim()) {
    flash(t('web.synergy.err_required'))
    return
  }
  if (drugA.value.trim().toLowerCase() === drugB.value.trim().toLowerCase()) {
    // 同药组合(2026-11-18 定案:仅警告不拦截):引擎照算(恒等式 s(A,A)=−d_AA/n
    // 有单测),输出为平凡重叠;真需要的用户(如自检基线)可跑,结果面板另有徽章
    flash(t('web.synergy.warn_same_drug'))
  }
  void job.run(api.synergyAnalyze({
    drug_a: drugA.value.trim(),
    drug_b: drugB.value.trim(),
    disease: disease.value.trim(),
    top_n: topN.value,
    n_perm: nPerm.value
  }))
}

/* 后端错误按稳定令牌映射 web.synergy.*(api/synergy.py 约定);未识别令牌原样透出 */
const errorMessage = computed(() => {
  const e = job.error.value
  if (!e) return ''
  const map: [string, string][] = [
    ['too_few_targets_a', t('web.synergy.err_too_few_a')],
    ['too_few_targets_b', t('web.synergy.err_too_few_b')],
    ['no_targets_a', t('web.synergy.err_no_targets_a')],
    ['no_targets_b', t('web.synergy.err_no_targets_b')],
    ['no_targets_disease', t('web.synergy.err_no_targets_disease')],
    ['InteractomeNotBuiltError', t('web.synergy.err_not_built')]
  ]
  for (const [token, text] of map) {
    if (e.includes(token)) return text
  }
  return e
})

/* 进度文案:优先 stage 码映射(核心只发 fraction,不带文案),兜底后端 text */
const stageText = computed(() => {
  const stage = String((job.extra.value as Record<string, unknown>).stage ?? '')
  if (stage === 'targets') return t('web.synergy.stage_targets')
  if (stage === 'separation') return t('web.synergy.stage_separation')
  if (stage === 'proximity') return t('web.synergy.stage_proximity')
  return job.text.value
})

/* 归档下载:走 downloads 集合的打包端点,拿到的就是该次运行目录的 ZIP */
async function downloadResult(): Promise<void> {
  if (!archiveDir.value) return
  const blob = await api.zip('downloads', [archiveDir.value])
  downloadBlob(blob, `${archiveDir.value}.zip`)
}
</script>

<template>
  <p class="kicker s1">{{ t('web.synergy.kicker') }}</p>
  <h1 class="title s2">{{ t('web.synergy.title') }}</h1>
  <p class="sub s3">{{ t('web.synergy.sub') }}</p>

  <div class="panel s4">
    <datalist id="synergy-drug-options">
      <option v-for="d in drugOptions" :key="d" :value="d"></option>
    </datalist>
    <datalist id="synergy-disease-options">
      <option v-for="d in diseaseOptions" :key="d" :value="d"></option>
    </datalist>

    <div class="names-form">
      <label class="field">
        <input v-model="drugA" spellcheck="false" list="synergy-drug-options"
               autocomplete="off" :placeholder="t('web.synergy.form_drug_a')" />
        <i class="line"></i>
      </label>
      <span class="cross">×</span>
      <label class="field">
        <input v-model="drugB" spellcheck="false" list="synergy-drug-options"
               autocomplete="off" :placeholder="t('web.synergy.form_drug_b')" />
        <i class="line"></i>
      </label>
    </div>

    <div class="disease-block">
      <span class="ctl-label">{{ t('web.synergy.form_disease') }}</span>
      <label class="field">
        <input v-model="disease" spellcheck="false" list="synergy-disease-options"
               autocomplete="off" :placeholder="t('web.synergy.ph_disease')" />
        <i class="line"></i>
      </label>
    </div>

    <div class="param-row">
      <span class="ctl-label">{{ t('web.synergy.param_n_label') }}</span>
      <span class="lv-group">
        <button v-for="opt in N_OPTIONS" :key="String(opt)" class="lv" type="button"
                :class="{ on: topN === opt }" @click="topN = opt">
          {{ opt === null ? t('web.synergy.param_n_full') : opt }}
        </button>
      </span>
      <span class="hint">{{ t('web.synergy.param_n_hint') }}</span>
    </div>
    <div class="param-row">
      <span class="ctl-label">{{ t('web.synergy.param_perm_label') }}</span>
      <span class="lv-group">
        <button v-for="opt in PERM_OPTIONS" :key="opt" class="lv" type="button"
                :class="{ on: nPerm === opt }" @click="nPerm = opt">{{ opt }}</button>
      </span>
      <span class="hint">{{ t('web.synergy.param_perm_hint') }}</span>
    </div>

    <div class="xlinks">
      <button class="bar-btn strong" type="button" :disabled="job.running.value" @click="run">
        {{ t('web.synergy.run_btn') }}
      </button>
      <a class="bar-btn" href="#/synergy-results">{{ t('web.nav.synergy_results_link') }}</a>
    </div>

    <div v-if="job.running.value" class="progress">
      <div class="progress-track">
        <div class="progress-fill"
             :style="{ width: `${((job.fraction.value ?? 0) * 100).toFixed(0)}%` }"></div>
      </div>
      <p class="progress-text">{{ stageText || t('web.synergy.stage_targets') }}</p>
    </div>
    <p v-if="errorMessage" class="notice err">{{ errorMessage }}</p>
    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>

  <div v-if="result" class="panel s5">
    <SynergyResultPanel :result="result" />

    <div v-if="archiveDir" class="xlinks arch">
      <button class="bar-btn" type="button" @click="downloadResult">
        {{ t('web.synergy.res_download') }}
      </button>
      <span class="hint">{{ t('web.synergy.res_archived') }}</span>
    </div>
  </div>
</template>

<style scoped>
.names-form {
  display: flex;
  align-items: flex-end;
  gap: 18px;
}
.names-form .field {
  flex: 1;
}
.disease-block {
  margin-top: 18px;
}
.disease-block .ctl-label {
  display: block;
  margin-bottom: 10px;
}
.param-row {
  display: flex;
  align-items: baseline;
  gap: 12px;
  margin-top: 16px;
  flex-wrap: wrap;
}
.hint {
  font-size: 10.5px;
  letter-spacing: 0.06em;
  color: var(--c-faint);
}
.xlinks {
  margin-top: 20px;
}
.xlinks.arch {
  margin-top: 18px;
}
</style>
