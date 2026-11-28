<script setup lang="ts">
/** 设置(App.py 齿轮职责的 TS 化,四栏侧边栏,排版对齐旧 Analysis settings):
 * ① 数据来源(抓取规模/STRING 版本) ② 分析算法(降噪阈值/ChEMBL 黑名单)
 * ③ 分子对接(口袋阈值/Vina·Open Babel 路径) ④ 通用设置(密钥/CLI 产物目录)。
 * 密钥与路径走 save_api_env(白名单+消毒+set_key+运行时回填);分析参数走
 * /api/settings/ui(后端值变化时自动清 merged 整合缓存并同步互作组置信度)。 */
import { onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { api, type SettingsSnapshot } from '../api'
import SettingsAlgoSection from '../components/SettingsAlgoSection.vue'
import SettingsDockingSection from '../components/SettingsDockingSection.vue'
import SettingsGeneralSection from '../components/SettingsGeneralSection.vue'
import SettingsSourcesSection from '../components/SettingsSourcesSection.vue'
import { useNotice } from '../composables/useNotice'

const { t } = useI18n()
const { notice, flash } = useNotice(3600)

/* ---------------- 侧边栏 ---------------- */

const SECTIONS = [
  { key: 'sources', label: 'web.settings.sec_sources' },
  { key: 'algo', label: 'web.settings.sec_algo' },
  { key: 'docking', label: 'web.settings.sec_docking' },
  { key: 'general', label: 'web.settings.sec_general' }
] as const
type SectionKey = (typeof SECTIONS)[number]['key']
const active = ref<SectionKey>('sources')

/* 快照单一真相源:四栏组件 watch(props.snapshot) 回填各自 ref;保存成功
 * emit('saved', 新快照) 也走这里统一更新,④ 栏无载荷(走 load 重新拉取,
 * 因 save_api_env 不返回快照)。栏内未保存的修改不丢——KeepAlive 缓存实例。 */
const snapshot = ref<SettingsSnapshot | null>(null)

async function load(): Promise<void> {
  snapshot.value = await api.settings()
}

function onSaved(s?: SettingsSnapshot): void {
  s ? (snapshot.value = s) : void load()
}

onMounted(() => void load().catch(() => flash(t('web.settings.read_fail'))))

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

        <!-- 四栏组件(KeepAlive 缓存实例:未保存的修改切栏不丢,同拆分前 ref 存活语义) -->
        <KeepAlive>
          <!-- ============ ① 数据来源 ============ -->
          <SettingsSourcesSection v-if="active === 'sources'" :snapshot="snapshot" @saved="onSaved" />
          <!-- ============ ② 分析算法 ============ -->
          <SettingsAlgoSection v-else-if="active === 'algo'" :snapshot="snapshot" @saved="onSaved" />
          <!-- ============ ③ 分子对接 ============ -->
          <SettingsDockingSection v-else-if="active === 'docking'" :snapshot="snapshot" @saved="onSaved" />
          <!-- ============ ④ 通用设置 ============ -->
          <SettingsGeneralSection v-else :snapshot="snapshot" :flash="flash" :reload="load" @saved="onSaved" />
        </KeepAlive>
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

/* 栏内容已拆入 Settings*Section 子组件:.sec/.opt-grid/.chip/.save-note 等
 * 基础控件类上提 screen.css(供四栏共用);.field/.lv/.bar-btn 相对全站基础版
 * 的紧凑变体(尺寸/间距)经 :deep 穿透到子组件,不影响其他页面的全局版。 */
.content :deep(.field) {
  width: 240px;
}
.content :deep(.field.wide) {
  width: min(420px, 100%);
}
.content :deep(.field input) {
  padding: 8px 2px;
  font-size: 13.5px;
}
.content :deep(.lv) {
  padding: 4px 0;
}
.content :deep(.bar-btn) {
  margin-top: 18px;
}
.content :deep(.add-row .bar-btn),
.content :deep(.ver-row + .set-row .bar-btn) {
  margin-top: 0;
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
