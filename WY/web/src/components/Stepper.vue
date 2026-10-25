<script setup lang="ts">
/** ± 步进器(齿轮同款交互):只允许按钮加减,到达上下限对应按钮置灰。 */
const props = withDefaults(defineProps<{
  modelValue: number
  min: number
  max: number
  step?: number
  /** 数值保留小数位(默认 2;富集 p 值等 3 位小数参数传 3) */
  digits?: number
}>(), { step: 1, digits: 2 })

const emit = defineEmits<{ (e: 'update:modelValue', v: number): void }>()

function dec(): void {
  emit('update:modelValue', Math.max(props.min, +(props.modelValue - props.step).toFixed(props.digits)))
}

function inc(): void {
  emit('update:modelValue', Math.min(props.max, +(props.modelValue + props.step).toFixed(props.digits)))
}
</script>

<template>
  <div class="stepper">
    <button class="st-btn" :disabled="modelValue <= min" @click="dec">−</button>
    <span class="st-val">{{ modelValue }}</span>
    <button class="st-btn" :disabled="modelValue >= max" @click="inc">+</button>
  </div>
</template>

<style scoped>
.stepper {
  display: inline-flex;
  align-items: center;
  gap: 12px;
}
.st-btn {
  width: 26px;
  height: 26px;
  border: 1px solid #e6e6e4;
  background: none;
  border-radius: 50%;
  font-size: 14px;
  color: #1a1a1a;
  cursor: pointer;
  transition: border-color 0.25s ease, opacity 0.25s ease;
}
.st-btn:hover:not(:disabled) {
  border-color: #1a1a1a;
}
.st-btn:disabled {
  opacity: 0.35;
  cursor: default;
}
.st-val {
  min-width: 34px;
  text-align: center;
  font-size: 14px;
  font-weight: 500;
  font-variant-numeric: tabular-nums;
  color: #1a1a1a;
}
</style>
