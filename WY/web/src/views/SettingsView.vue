<script setup lang="ts">
/** 设置(App.py 齿轮职责的 TS 化,四栏侧边栏,排版对齐旧 Analysis settings):
 * ① 数据来源(抓取规模/STRING 版本) ② 分析算法(降噪阈值/ChEMBL 黑名单)
 * ③ 分子对接(口袋阈值/Vina·Open Babel 路径) ④ 通用设置(密钥/CLI 产物目录)。
 * 密钥与路径走 save_api_env(白名单+消毒+set_key+运行时回填);分析参数走
 * /api/settings/ui(后端值变化时自动清 merged 整合缓存并同步互作组置信度)。 */
import { computed, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api } from '../api'
import Stepper from '../components/Stepper.vue'
import { useNotice } from '../composables/useNotice'
import { applyScale, persistScale, snapScale } from '../utils/fontScale'
import { applyTheme, normalizeTheme, persistTheme, type Theme } from '../utils/theme'

const { t } = useI18n()
const { notice, flash } = useNotice(3600)

type Snapshot = Awaited<ReturnType<typeof api.settings>>

/* ---------------- 侧边栏 ---------------- */

const SECTIONS = [
  { key: 'sources', label: 'web.settings.sec_sources' },
  { key: 'algo', label: 'web.settings.sec_algo' },
  { key: 'docking', label: 'web.settings.sec_docking' },
  { key: 'general', label: 'web.settings.sec_general' }
] as const
type SectionKey = (typeof SECTIONS)[number]['key']
const active = ref<SectionKey>('sources')

/* ---------------- 状态(一次快照全部载入) ---------------- */

const keys = ref({ disgenet: false, ncbi: false })
const dgKey = ref('')
const ncbiKey = ref('')
const cliResultsDir = ref('')

/* 字体大小(界面等比缩放):滑动条 80%–150%、步长 10%,savedScale 用于保存按钮的未修改判定 */
const FONT_SCALE_STEP = 0.1
const fontScale = ref(1)
const savedScale = ref(1)
const fontScaleDirty = computed(() => fontScale.value !== savedScale.value)

/* 颜色变换:浅色/深色(样品「墨黑」),savedTheme 用于保存按钮的未修改判定 */
const themeSel = ref<Theme>('light')
const savedTheme = ref<Theme>('light')
const themeDirty = computed(() => themeSel.value !== savedTheme.value)

const otTopN = ref(500)
const disgenetPages = ref(5)
const stringVersion = ref('')
const knownVersions = ref<string[]>([])
const customMode = ref(false)
const customVersion = ref('')

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

const pocketCutoff = ref(4.0)
const vinaDir = ref('')
const obabelDir = ref('')

function applySnapshot(s: Snapshot): void {
  keys.value = s.keys
  cliResultsDir.value = s.cli_results_dir
  fontScale.value = snapScale(s.ui_font_scale)
  savedScale.value = fontScale.value
  themeSel.value = normalizeTheme(s.ui_theme)
  savedTheme.value = themeSel.value
  otTopN.value = s.ot_top_n
  disgenetPages.value = s.disgenet_pages
  stringVersion.value = s.string_version
  knownVersions.value = s.known_versions
  blacklist.value = s.chembl_blacklist
  blacklistOptions.value = s.chembl_blacklist_options
  pchembl.value = s.chembl_pchembl_threshold
  dgidbMin.value = s.dgidb_min_score
  corroboration.value = s.chembl_corroboration_required
  otMin.value = s.ot_min_score
  enrichP.value = s.enrichment_pvalue
  louvain.value = s.louvain_resolution
  ppiConf.value = s.ppi_confidence
  pocketCutoff.value = s.pocket_cutoff
  vinaDir.value = s.vina_dir
  obabelDir.value = s.obabel_dir
}

async function load(): Promise<void> {
  applySnapshot(await api.settings())
}

onMounted(() => void load().catch(() => flash(t('web.settings.read_fail'))))

/* 保存后的统一提示:基础文案 + 整合缓存清理说明 */
function flashSaved(purged: number): void {
  flash(purged > 0
    ? `${t('web.settings.saved')} ${t('web.settings.purged', { n: purged })}`
    : t('web.settings.saved'))
}

/* ---------------- ① 数据来源 ---------------- */

async function saveSources(): Promise<void> {
  try {
    const r = await api.saveUiSettings({
      OT_TOP_N: otTopN.value, DISGENET_PAGES: disgenetPages.value
    })
    applySnapshot(r.snapshot)
    flashSaved(r.purged_merged)
  } catch (e) {
    flash(e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}

async function applyStringVersion(): Promise<void> {
  const target = customMode.value ? customVersion.value.trim() : stringVersion.value
  if (!target) return
  try {
    const r = await api.setStringVersion(target)
    flash(r.changed ? t('web.settings.switched', { v: r.version }) : t('web.settings.no_change'))
    stringVersion.value = r.version
  } catch (e) {
    flash(e instanceof Error ? e.message : t('web.settings.switch_fail'))
  }
}

function chooseVersion(v: string): void {
  stringVersion.value = v
}

/* ---------------- ② 分析算法 ---------------- */

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
  try {
    const r = await api.saveUiSettings({
      CHEMBL_TARGET_BLACKLIST: blacklist.value,
      CHEMBL_PCHEMBL_THRESHOLD: pchembl.value,
      DGIDB_MIN_SCORE: dgidbMin.value,
      CHEMBL_CORROBORATION_REQUIRED: corroboration.value,
      OT_MIN_SCORE: otMin.value,
      ENRICHMENT_PVALUE: enrichP.value,
      LOUVAIN_RESOLUTION: louvain.value,
      PPI_CONFIDENCE: ppiConf.value
    })
    applySnapshot(r.snapshot)
    flashSaved(r.purged_merged)
  } catch (e) {
    flash(e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}

/* ---------------- ③ 分子对接 ---------------- */

async function saveDocking(): Promise<void> {
  try {
    const ui = await api.saveUiSettings({ POCKET_CUTOFF: pocketCutoff.value })
    applySnapshot(ui.snapshot)
    const applied = (await api.saveApiEnv({
      VINA_DIR: vinaDir.value, OBABEL_DIR: obabelDir.value
    })).applied
    flashSaved(0)
    if (applied.length) flash(`${t('web.settings.saved')} · ${applied.join(', ')}`)
  } catch (e) {
    flash(e instanceof Error ? e.message : t('web.settings.save_fail'))
  }
}

/* ---------------- ④ 通用设置 ---------------- */

async function saveGeneral(): Promise<void> {
  try {
    if (fontScaleDirty.value || themeDirty.value) {
      // 缩放/主题即时生效 + 本地镜像,再随快照确认后端值(纯显示偏好,不触发缓存清理)
      const ui = await api.saveUiSettings({
        ...(fontScaleDirty.value ? { UI_FONT_SCALE: fontScale.value } : {}),
        ...(themeDirty.value ? { UI_THEME: themeSel.value } : {})
      })
      applySnapshot(ui.snapshot)
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
    flash(applied.length ? t('web.settings.keys_saved', { names: applied.join(', ') })
                         : t('web.settings.saved'))
    await load()
  } catch {
    flash(t('web.settings.keys_fail'))
  }
}
</script>

<template>
  <p class="kicker s1">{{ t('web.settings.kicker') }}</p>
  <h1 class="title s2">{{ t('web.settings.title') }}</h1>
  <p class="sub s3">{{ t('web.settings.sub') }}</p>

  <div class="panel wide s4">
    <p class="hint-line top-hint">{{ t('web.settings.page_hint') }}</p>

    <div class="layout">
      <!-- 侧边栏 -->
      <aside class="side">
        <button v-for="s in SECTIONS" :key="s.key" type="button" class="side-item"
                :class="{ on: active === s.key }" @click="active = s.key">
          <i class="side-dot" aria-hidden="true"></i>{{ t(s.label) }}
        </button>
      </aside>

      <!-- 内容区 -->
      <div class="content">
        <Transition name="note">
          <p v-if="notice" class="notice">{{ notice }}</p>
        </Transition>

        <!-- ============ ① 数据来源 ============ -->
        <template v-if="active === 'sources'">
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
        </template>

        <!-- ============ ② 分析算法 ============ -->
        <template v-else-if="active === 'algo'">
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

          <button class="bar-btn" @click="saveAlgo">{{ t('web.settings.save') }}</button>
        </template>

        <!-- ============ ③ 分子对接 ============ -->
        <template v-else-if="active === 'docking'">
          <div class="opt-grid">
            <label class="opt">
              <span class="opt-label">{{ t('web.settings.pocket_label') }}</span>
              <Stepper v-model="pocketCutoff" :min="3" :max="6" :step="0.5" />
            </label>
          </div>
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
        </template>

        <!-- ============ ④ 通用设置 ============ -->
        <template v-else>
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
        </template>
      </div>
    </div>
  </div>
</template>

<style scoped>
.top-hint {
  margin: 0 0 22px;
}
.layout {
  display: flex;
  align-items: flex-start;
  gap: 48px;
}
/* 侧边栏:旧版单选列表的极简化(圆点指示,文字开关) */
.side {
  flex: none;
  width: 132px;
  display: flex;
  flex-direction: column;
  gap: 2px;
  position: sticky;
  top: 0;
}
.side-item {
  display: flex;
  align-items: center;
  gap: 9px;
  border: none;
  background: none;
  padding: 7px 0;
  font-family: inherit;
  font-size: 12.5px;
  letter-spacing: 0.1em;
  color: var(--c-faint);
  cursor: pointer;
  transition: color 0.25s ease;
  text-align: left;
}
.side-item:hover {
  color: var(--c-ink2);
}
.side-item.on {
  color: var(--c-ink);
  font-weight: 500;
}
.side-dot {
  flex: none;
  width: 5px;
  height: 5px;
  border-radius: 50%;
  border: 1px solid var(--c-dot);
  transition: background 0.25s ease, border-color 0.25s ease;
}
.side-item.on .side-dot {
  background: var(--c-accent);
  border-color: var(--c-ink);
}
.content {
  flex: 1;
  min-width: 0;
}

/* 参数网格:标签在上、控件在下(AnalysisView 同款排版) */
.opt-grid {
  display: flex;
  flex-wrap: wrap;
  gap: 22px 44px;
}
.opt {
  display: flex;
  flex-direction: column;
  gap: 6px;
  min-width: 150px;
}
.opt-label {
  font-size: 10.5px;
  letter-spacing: 0.14em;
  color: var(--c-ink3);
}
.corro-row {
  margin-top: 16px;
}

/* 黑名单标签 */
.chips {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin-top: 10px;
}
.chip {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  border: 1px solid var(--c-line);
  padding: 3px 9px;
  font-size: 11.5px;
  letter-spacing: 0.06em;
  color: var(--c-ink);
  font-variant-numeric: tabular-nums;
}
button.chip {
  cursor: pointer;
  color: var(--c-ink2);
  transition: color 0.25s ease, border-color 0.25s ease;
}
button.chip:hover {
  color: var(--c-ink);
  border-color: var(--c-ink);
}
.chip-x {
  border: none;
  background: none;
  padding: 0;
  font-size: 13px;
  line-height: 1;
  color: var(--c-sub);
  cursor: pointer;
  transition: color 0.25s ease;
}
.chip-x:hover {
  color: var(--c-danger);
}
.add-row {
  display: flex;
  align-items: baseline;
  gap: 18px;
  margin-top: 14px;
}
.add-row .field {
  width: 200px;
}

/* 本节专用控件(其余复用全局 .sec/.set-row/.field/.lv/.bar-btn/.chk) */
.field {
  position: relative;
  width: 240px;
}
.field.wide {
  width: min(420px, 100%);
}
.field input {
  width: 100%;
  box-sizing: border-box;
  border: none;
  outline: none;
  background: transparent;
  border-bottom: 1px solid var(--c-line);
  padding: 8px 2px;
  font-size: 13.5px;
  font-family: inherit;
  color: var(--c-ink);
  border-radius: 0;
}
.line {
  position: absolute;
  left: 0;
  right: 0;
  bottom: 0;
  height: 1px;
  background: var(--c-accent);
  transform: scaleX(0);
  transition: transform 0.45s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.field:focus-within .line {
  transform: scaleX(1);
}
.state {
  font-size: 11px;
  letter-spacing: 0.08em;
  color: var(--c-ink4);
}
.state.on {
  color: var(--c-ink);
}
.ver-row {
  display: flex;
  align-items: baseline;
  gap: 12px;
  padding: 6px 0;
  flex-wrap: wrap;
}
/* 字体大小滑动条:发丝线轨道 + 圆形黑滑钮,右侧等宽百分比读数 */
.range-row {
  display: flex;
  align-items: center;
  gap: 16px;
  padding: 10px 0 2px;
  max-width: 380px;
}
.range-row input[type="range"] {
  flex: 1;
  appearance: none;
  -webkit-appearance: none;
  height: 2px;
  background: var(--c-line);
  outline: none;
  cursor: pointer;
}
.range-row input[type="range"]::-webkit-slider-thumb {
  -webkit-appearance: none;
  appearance: none;
  width: 14px;
  height: 14px;
  border-radius: 50%;
  background: var(--c-accent);
  border: none;
}
.range-row input[type="range"]::-moz-range-thumb {
  width: 14px;
  height: 14px;
  border-radius: 50%;
  background: var(--c-accent);
  border: none;
}
.range-row input[type="range"]::-moz-range-track {
  height: 2px;
  background: var(--c-line);
}
.range-val {
  flex: none;
  min-width: 44px;
  text-align: right;
  font-size: 12.5px;
  letter-spacing: 0.06em;
  color: var(--c-ink);
  font-variant-numeric: tabular-nums;
}
.sec {
  margin-top: 26px;
  padding-top: 18px;
  border-top: 1px solid var(--c-line-soft);
}
.sec-title {
  margin: 0 0 8px;
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0.14em;
  color: var(--c-ink);
}
.bar-btn {
  margin-top: 18px;
  border: none;
  background: none;
  font-family: inherit;
  font-size: 12px;
  letter-spacing: 0.08em;
  color: var(--c-ink);
  cursor: pointer;
  padding: 4px 0;
  border-bottom: 1px solid var(--c-accent);
  transition: opacity 0.25s ease;
}
.add-row .bar-btn,
.ver-row + .set-row .bar-btn {
  margin-top: 0;
}
.bar-btn:disabled {
  cursor: default;
  border-bottom-color: var(--c-ghost);
  color: var(--c-ink4);
}
.lv {
  border: none;
  background: none;
  padding: 4px 0;
  font-family: inherit;
  font-size: 12.5px;
  letter-spacing: 0.12em;
  color: var(--c-faint);
  cursor: pointer;
  transition: color 0.25s ease;
}
.lv:hover {
  color: var(--c-ink2);
}
.lv.on {
  color: var(--c-ink);
  font-weight: 500;
}
.notice {
  margin: 0 0 14px;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: var(--c-dim);
}

@media (max-width: 760px) {
  .layout {
    flex-direction: column;
    gap: 18px;
  }
  .side {
    width: 100%;
    flex-direction: row;
    flex-wrap: wrap;
    gap: 4px 22px;
  }
}
</style>
