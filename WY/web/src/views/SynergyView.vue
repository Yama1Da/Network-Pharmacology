<script setup lang="ts">
/** 协同预测(二期占位页):导航/路由先行,引擎落地后在此补药物对表单、
 * z_A/z_B 网络分离度计算、协同/拮抗判定与靶点浓度可达性标注列。
 * 深色适配:全页只用家族类(kicker/title/sub/panel/notice),色值全走
 * var(--c-*),跟随 :root[data-theme="dark"] 变量表,无自带 scoped 样式。 */
import { useI18n } from 'vue-i18n'

const { t } = useI18n()
</script>

<template>
  <p class="kicker s1">{{ t('web.synergy.kicker') }}</p>
  <h1 class="title s2">{{ t('web.synergy.title') }}</h1>
  <p class="sub s3">{{ t('web.synergy.sub') }}</p>

  <div class="panel s4">
    <p class="notice">{{ t('web.synergy.building') }}</p>
  </div>
</template>
