简体中文 | [English](README_EN.md)

# 网络药理学分析平台（Network Pharmacology Platform）

基于 Streamlit 的一站式网络药理学分析工具：输入药物与疾病名称，自动完成
「靶点获取 → 靶点整合 → PPI 网络构建 → 富集分析 → 分子对接」的完整流程。


## 功能特性

- **疾病/药物靶点获取**：接入 DisGeNET、DGIdb、Open Targets、ChEMBL 四个数据源，
  自动整合去重，并过滤非特异性结合蛋白（白蛋白、CYP 代谢酶等）降噪
- **PPI 网络**：基于 STRING v12.0 构建蛋白互作网络，支持多种可视化风格展示
- **富集分析**：GO / KEGG 等通路富集，g:Profiler 失败时自动切换 Enrichr
- **分子对接**：集成 AutoDock Vina，RDKit + Meeko 自动完成配体准备，内置 3D 结果查看 + 聚焦口袋（虚化周边）
- **工程化**：请求级缓存与 DisGeNET 本地限流（学术配额内）、中英双语界面、
  每日滚动日志、靶点/疾病数据管理页面

## 环境要求

- Python ≥ 3.11
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)（加入 PATH，或在 `WY/API.env` （请将API.env.example重命名为API.env）中指定 `VINA_DIR`）
- [DisGeNET](https://www.disgenet.org/) 学术账号 API Key（免费申请）/ 其他级别账号 API Key

## 快速开始

```bash
# 1. 创建虚拟环境
python -m venv .venv
.venv\Scripts\activate        # Windows（Linux/macOS: source .venv/bin/activate）

# 2. 安装依赖
pip install -r requirements.txt
建议 conda install -c conda-forge openbabel #`openbabel>=3.2.0` 可以 pip 安装，但 conda 安装更省事

# 3. 配置 API.env
编辑 API.env，填入 DisGeNET_API_KEY 与 VINA_DIR 路径

# 4. 启动
streamlit run App.py  #终端输入
```

## 使用流程

1. **首页发起分析**：输入药物名与疾病名（默认示例为 metformin 与 Type 2 diabetes mellitus）,现暂时不支持中文疾病名输入，
   点击「开始分析」。程序会并发请求 ChEMBL / DGIdb / Open Targets / DisGeNET 四个数据源。
2. **查看靶点概览**：分析完成后首页展示药物靶点数、疾病靶点数、交集靶点数与韦恩图，
   展开折叠面板可查看各数据源的靶点明细表（按来源数与打分排序）。
3. **数据管理**：通过首页左上角的「数据管理」菜单进入各信息页面，
   查看、增删或修正抓取到的疾病靶点与药物靶点。
4. **PPI 网络与富集**：从菜单进入「数据网络」页面，基于整合后的靶点构建 STRING
   互作网络（可选多种可视化风格），并执行 GO / KEGG 富集分析。
5. **分子对接**：从菜单进入「分子对接」页面，选择核心靶点进行 AutoDock Vina 对接，
   在线查看 3D 结合构象与聚焦口袋效果。
6. **日志与设置**：右上角可切换中/英文界面并查看运行日志；左下角齿轮可调整全局设置
   （如靶点过滤阈值），设置变更后整合缓存自动失效重建。

所有中间结果自动缓存（默认 6 小时，见 `WY/config.py`），重复分析秒级返回。

## 配置说明

| 文件 | 用途 |
|---|---|
| `WY/API.env` | 密钥与可执行文件路径（`DisGeNET_API_KEY`、`VINA_DIR`），不入库 |
| `WY/config.py` | 数据源地址、限流/重试、缓存 TTL、靶点过滤阈值 |

> 改行为看 `WY/config.py`，填密钥看 `WY/API.env`：`WY/API.env` 只放不提交到 git 的敏感信息
> 与本机路径，其余请求控制、缓存有效期、靶点降噪过滤阈值集中在 `WY/config.py`。

## 常见问题

- **DisGeNET 返回 429？** 程序已在学术配额（50 次/分）内做了本地限流，仍触发时稍等即可；
  配额为每天 25,000 次，正常单次分析仅消耗 1 次请求。如有其他等级配额，请自行更改 `WY/config.py`中的
  `DISGENET_RATE_LIMIT_PER_MINUTE` `DISGENET_RATE_LIMIT_PER_DAY` `DISGENET_RATE_LIMIT_PER_MONTH`。
- **富集分析没结果？** g:Profiler 服务不可达时会自动熔断 30 分钟并改用 Enrichr。
- **改了过滤阈值后结果没变？** 靶点过滤阈值变更需同步递增 `TARGET_FILTER_VERSION`，
  旧整合缓存会自动失效重建。
- **分子对接找不到 Vina？** 确认 `vina` 已加入 PATH 或环境变量存在 `VINA_DIR` ，或在 `WY/API.env` 中将 `VINA_DIR`
  指向 Vina 所在目录（也可填可执行文件全路径）。

## 免责声明

本项目仅供科研用途，分析结果请以官方数据库与实验验证为准。
使用数据时请遵循各数据源（STRING、DisGeNET、ChEMBL、DGIdb、Open Targets）的引用要求。