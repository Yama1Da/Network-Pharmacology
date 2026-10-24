import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

// 开发期:vite(5175)只管前端热更新,/api 代理给 FastAPI(8000,先 npap serve --no-browser);
// 生产(npm run build)产出 dist/ 由 FastAPI 同源伺服,此后与 Vite、Node 再无关系。
export default defineConfig({
  plugins: [vue()],
  server: {
    port: 5175,
    host: true,
    proxy: {
      '/api': 'http://127.0.0.1:8000'
    }
  }
})
