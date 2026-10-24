<script setup lang="ts">
/** 设置(App.py 齿轮职责的 TS 化,阶段 3 补全):
 * 密钥 write-only / 工具路径明文 / STRING 版本(探活+确认) / 口袋阈值。
 * 密钥与路径走 core.management.save_api_env(白名单+消毒+set_key+运行时回填);
 * 其余齿轮参数(抓取规模/置信度等)暂留 Streamlit 齿轮,换届逐步开放。 */
import { onMounted, ref } from 'vue'
import { api } from '../api'
import Stepper from '../components/Stepper.vue'

const keys = ref({ disgenet: false, ncbi: false })
const dgKey = ref('')
const ncbiKey = ref('')
const vinaDir = ref('')
const obabelDir = ref('')
const cliResultsDir = ref('')
const stringVersion = ref('')
const knownVersions = ref<string[]>([])
const customMode = ref(false)
const customVersion = ref('')
const pocketCutoff = ref(4.0)
const notice = ref('')
let noticeTimer = 0

function flash(msg: string): void {
  notice.value = msg
  window.clearTimeout(noticeTimer)
  noticeTimer = window.setTimeout(() => (notice.value = ''), 3200)
}

async function load(): Promise<void> {
  const s = await api.settings()
  keys.value = s.keys
  vinaDir.value = s.vina_dir
  obabelDir.value = s.obabel_dir
  cliResultsDir.value = s.cli_results_dir
  stringVersion.value = s.string_version
  knownVersions.value = s.known_versions
  pocketCutoff.value = s.pocket_cutoff
}

onMounted(() => void load().catch(() => flash('设置读取失败')))

async function saveKeys(): Promise<void> {
  const changes: Record<string, string> = {}
  if (dgKey.value.trim()) changes.DisGeNET_API_KEY = dgKey.value
  if (ncbiKey.value.trim()) changes.NCBI_API_KEY = ncbiKey.value
  if (!Object.keys(changes).length) return
  try {
    const applied = (await api.saveApiEnv(changes)).applied
    flash(`已写入 API.env:${applied.join(', ')}(即时生效)`)
    dgKey.value = ''
    ncbiKey.value = ''
    await load()
  } catch {
    flash('写入失败')
  }
}

async function savePaths(): Promise<void> {
  try {
    const applied = (await api.saveApiEnv({
      VINA_DIR: vinaDir.value,
      OBABEL_DIR: obabelDir.value,
      CLI_RESULTS_DIR: cliResultsDir.value
    })).applied
    flash(applied.length ? `已生效:${applied.join(', ')}(可执行文件已重探)` : '无变更')
    await load()
  } catch {
    flash('保存失败')
  }
}

async function applyStringVersion(): Promise<void> {
  const target = customMode.value ? customVersion.value.trim() : stringVersion.value
  if (!target) return
  try {
    const r = await api.setStringVersion(target)
    flash(r.changed ? `已切换到 ${r.version}(互作组缓存按版本隔离,已清理)` : '版本未变化')
    stringVersion.value = r.version
  } catch (e) {
    flash(e instanceof Error ? e.message : '切换失败')
  }
}

async function savePocket(): Promise<void> {
  try {
    const r = await api.saveUiSettings({ POCKET_CUTOFF: pocketCutoff.value })
    pocketCutoff.value = r.pocket_cutoff
    flash('口袋判定阈值已保存并即时生效')
  } catch {
    flash('保存失败')
  }
}

function chooseVersion(v: string): void {
  stringVersion.value = v
}
</script>

<template>
  <p class="kicker s1">SETTINGS</p>
  <h1 class="title s2">设置</h1>
  <p class="sub s3">密钥与路径写入 API.env&nbsp;·&nbsp;改动即时生效</p>

  <div class="panel s4">
    <!-- API 密钥(write-only) -->
    <section class="sec">
      <h3 class="sec-title">API 密钥</h3>
      <p class="hint-line">不回显当前值;输入新值即覆盖,留空 = 不修改。</p>
      <div class="set-row">
        <label class="set-label">DisGeNET API Key</label>
        <div class="field"><input v-model="dgKey" type="password" placeholder="********"
                                  spellcheck="false" /><i class="line"></i></div>
        <span class="state" :class="{ on: keys.disgenet }">{{ keys.disgenet ? '已配置' : '未配置' }}</span>
      </div>
      <div class="set-row">
        <label class="set-label">NCBI API Key</label>
        <div class="field"><input v-model="ncbiKey" type="password" placeholder="********"
                                  spellcheck="false" /><i class="line"></i></div>
        <span class="state" :class="{ on: keys.ncbi }">{{ keys.ncbi ? '已配置' : '未配置' }}</span>
      </div>
      <button class="bar-btn" :disabled="!dgKey.trim() && !ncbiKey.trim()" @click="saveKeys">保存密钥</button>
    </section>

    <!-- 工具路径(明文) -->
    <section class="sec">
      <h3 class="sec-title">本机工具路径</h3>
      <p class="hint-line">路径非密钥,可明文查看;留空 = 不修改(回退自动探测请手动编辑 API.env 后重启)。</p>
      <div class="set-row">
        <label class="set-label">Vina 目录</label>
        <div class="field wide"><input v-model="vinaDir" spellcheck="false" /><i class="line"></i></div>
      </div>
      <div class="set-row">
        <label class="set-label">Open Babel 目录</label>
        <div class="field wide"><input v-model="obabelDir" spellcheck="false" /><i class="line"></i></div>
      </div>
      <div class="set-row">
        <label class="set-label">CLI 结果目录</label>
        <div class="field wide"><input v-model="cliResultsDir" spellcheck="false" /><i class="line"></i></div>
      </div>
      <button class="bar-btn" @click="savePaths">保存路径</button>
    </section>

    <!-- STRING 版本 -->
    <section class="sec">
      <h3 class="sec-title">STRING 数据版本(当前 {{ stringVersion }})</h3>
      <p class="hint-line">查询 API 与背景互作组按版本同源派生;切换前自动探测下载站可用性。</p>
      <div class="ver-row">
        <button v-for="v in knownVersions" :key="v" class="lv" :class="{ on: stringVersion === v }"
                @click="chooseVersion(v)">{{ v }}</button>
        <button class="lv" :class="{ on: customMode }" @click="customMode = !customMode">自定义…</button>
      </div>
      <div v-if="customMode" class="set-row">
        <label class="set-label">版本号</label>
        <div class="field"><input v-model="customVersion" placeholder="v13.0" spellcheck="false" /><i class="line"></i></div>
        <button class="bar-btn" :disabled="!customVersion.trim()" @click="applyStringVersion">切换到该版本</button>
      </div>
      <div v-else-if="stringVersion && !knownVersions.includes(stringVersion)" class="ver-row">
        <span class="hint-line">当前版本 {{ stringVersion }} 不在已知清单:重新选择上面某版即可切换回去。</span>
      </div>
      <button v-else class="bar-btn" :disabled="true">当前版本 · 切换请选择上方其他版本</button>
    </section>

    <!-- 口袋阈值 -->
    <section class="sec">
      <h3 class="sec-title">分子对接</h3>
      <div class="set-row">
        <label class="set-label">口袋判定阈值(Å)</label>
        <Stepper v-model="pocketCutoff" :min="3" :max="6" :step="0.5" />
        <button class="bar-btn" @click="savePocket">保存</button>
      </div>
    </section>

    <Transition name="note">
      <p v-if="notice" class="notice">{{ notice }}</p>
    </Transition>
  </div>
</template>

<style scoped>
.set-row {
  display: flex;
  align-items: baseline;
  gap: 18px;
  padding: 10px 0;
}
.set-label {
  flex: none;
  width: 130px;
  font-size: 11.5px;
  letter-spacing: 0.1em;
  color: #8a8a88;
}
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
  border-bottom: 1px solid #e6e6e4;
  padding: 8px 2px;
  font-size: 13.5px;
  font-family: inherit;
  color: #1a1a1a;
  border-radius: 0;
}
.line {
  position: absolute;
  left: 0;
  right: 0;
  bottom: 0;
  height: 1px;
  background: #1a1a1a;
  transform: scaleX(0);
  transition: transform 0.45s cubic-bezier(0.22, 0.61, 0.36, 1);
}
.field:focus-within .line {
  transform: scaleX(1);
}
.state {
  font-size: 11px;
  letter-spacing: 0.08em;
  color: #c6c6c4;
}
.state.on {
  color: #1a1a1a;
}
.ver-row {
  display: flex;
  align-items: baseline;
  gap: 12px;
  padding: 6px 0;
  flex-wrap: wrap;
}
.sec {
  margin-top: 26px;
  padding-top: 18px;
  border-top: 1px solid #f0f0ee;
}
.sec:first-of-type {
  border-top: none;
  margin-top: 0;
  padding-top: 0;
}
.sec-title {
  margin: 0 0 8px;
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0.14em;
  color: #1a1a1a;
}
.bar-btn {
  border: none;
  background: none;
  font-family: inherit;
  font-size: 12px;
  letter-spacing: 0.08em;
  color: #1a1a1a;
  cursor: pointer;
  padding: 4px 0;
  border-bottom: 1px solid #1a1a1a;
  transition: opacity 0.25s ease;
}
.bar-btn:disabled {
  cursor: default;
  border-bottom-color: #dedede;
  color: #c6c6c4;
}
.lv {
  border: none;
  background: none;
  padding: 4px 0;
  font-family: inherit;
  font-size: 12.5px;
  letter-spacing: 0.12em;
  color: #cfcfcd;
  cursor: pointer;
  transition: color 0.25s ease;
}
.lv:hover {
  color: #8a8a88;
}
.lv.on {
  color: #1a1a1a;
  font-weight: 500;
}
.notice {
  margin: 12px 0 0;
  font-size: 11px;
  letter-spacing: 0.08em;
  color: #9a9a98;
}
</style>
