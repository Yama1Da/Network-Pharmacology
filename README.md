简体中文 | [EN](README_EN.md)

# 网络药理学分析平台（NPAP）

[![CI](https://github.com/Yama1Da/Network-Pharmacology/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/Yama1Da/Network-Pharmacology/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.13%2C%20%3C3.14-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

一站式网络药理学分析平台：输入药物与疾病名称，自动完成
「靶点获取 → 靶点整合 → PPI 网络构建 → 富集分析 → 分子对接」完整流程，
并可结合 GEO 表达谱与多算法机器学习共识，收缩值得湿实验验证的靶点清单。
本项目主服务于网络药理学科研人员与本科生。
网页端全流程可视化，CLI 支持批量出表，不包含出图片。

## 功能特性

- **疾病/药物靶点获取**：接入 DisGeNET、DGIdb、Open Targets、ChEMBL 四个数据源，
  自动整合去重，并过滤非特异性结合蛋白（白蛋白、CYP 代谢酶等）降噪，过滤条件可调。
- **PPI 网络**：基于 STRING 构建蛋白互作网络（默认 v12.0，设置中可切换数据版本），
  拓扑分析给出 MCC 核心靶点，支持多种可视化风格。
- **靶点优先级**：本地缓存 STRING 全网络作背景互作组，支撑邻近度 Z-score 与
  RWR 全网络排序；数据管理页一键构建，镜像/官方双源自动回退 + 分段并发下载。
- **富集分析**：GO / KEGG 通路富集，g:Profiler 失败自动切换 Enrichr；双源均带
  熔断降级与自动恢复，都不可用时页面提示"约 X 分钟后自动重试"。
- **分子对接**：集成 AutoDock Vina，RDKit + Meeko 自动完成配体准备，
  内置 3D 结果查看与聚焦口袋（虚化周边）展示。
- **表达谱与 ML 共识筛选**：检索并下载 GEO 疾病表达数据集，自动完成差异表达
  分析（DEG）后，用 LASSO / Elastic Net / Ridge / SVM-RFE / 随机森林+Boruta /
  XGBoost 六算法交叉筛选，共识打分 + CV-ROC 评估，再与平台排序靶点取交集，
  输出湿实验优先级清单（shortlist）。
- **工程化**：请求级缓存、DisGeNET 本地限流（学术配额内）、富集双源熔断降级、
  中英双语界面、每日滚动日志、靶点数据管理页面。

## 环境要求

Python 依赖无需手动安装：启动脚本（Windows 双击 `start_npap.bat`，Linux/macOS
执行 `./start_npap.sh`）首次运行会自动创建虚拟环境并安装。最低要求：

- Python 3.13（精确版本范围见顶部 Python 徽章）
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/releases) 与
  [Open Babel](https://openbabel.org/)：仅「分子对接」功能需要。装完在 PATH 即被
  自动发现，也可在设置页或 `WY/API.env` 中指定 `VINA_DIR` / `OBABEL_DIR`
- [DisGeNET](https://www.disgenet.org/) API Key：疾病–基因关联数据库，疾病靶点
  来源之一，学术账号免费申请。在网页「设置」页填入即可，保存后自动写入
  `WY/API.env`；不填也能用本项目，仅 DisGeNET 靶点表为空，其余数据源不受影响

## 快速开始

第 1、2 步是必需的，完成后即可在网页跑通主流程；第 3、4 步按需。

### 1. 安装 Python 3.13

Windows 从 [python.org](https://www.python.org/downloads/) 安装，勾选
"Add python.exe to PATH"；macOS / Linux：

```bash
brew install python@3.13                       # macOS
sudo apt install python3.13 python3.13-venv    # Debian/Ubuntu
sudo dnf install python3.13                    # Fedora
```

### 2. 一键启动

Windows 双击 `start_npap.bat`（停止双击 `stop_npap.bat`）；
Linux/macOS 执行：
```bash
./start_npap.sh    # 停止: ./stop_npap.sh
```

首次运行自动创建虚拟环境并安装依赖（需几分钟，默认走清华 PyPI 镜像，
置空 `NPAP_PIP_INDEX` 环境变量可改用官方源），完成后自动打开浏览器
`http://127.0.0.1:8000`，没有打开就手动访问。

### 3. 安装分子对接工具（仅「分子对接」功能需要）

Windows 按「环境要求」里的链接下载安装包；
Linux/macOS：

```bash
conda install -c conda-forge vina openbabel
```

### 4. （可选）命令行批量分析

启动脚本已完成可编辑安装，激活虚拟环境后 `npap` 即全局可用：

```bash
source .venv/bin/activate   # Windows: .venv\Scripts\activate
npap --help                 # 详见下方「命令行批量分析（CLI）」一节
```

## 使用流程

1. **首页发起分析**：输入药物名与疾病名（默认示例 metformin 与
   Type 2 diabetes mellitus，暂仅支持英文疾病名），点击「开始分析」，
   程序并发请求四个数据源。
2. **查看靶点概览**：分析完成后展示药物/疾病/交集靶点数与韦恩图，
   展开折叠面板可查看各数据源靶点明细表（按来源数与打分排序）。
3. **数据管理**：左上角「数据管理」菜单进入，查看、增删或修正抓取到的
   疾病/药物靶点；同一页面可一键构建背景互作组（邻近度 / RWR 分析的前置）。
4. **数据网络**：基于整合后的靶点构建 STRING 互作网络（可选多种可视化风格），
   并执行 GO / KEGG 富集分析。
5. **分子对接**：选择核心靶点进行 AutoDock Vina 对接，在线查看 3D 结合
   构象与聚焦口袋效果。
6. **表达谱与 ML 共识筛选**（可选）：检索并选定 GEO 数据集（或上传本地
   表达矩阵），自动完成 DEG 与六算法共识筛选，下载 shortlist 作为湿实验
   优先级参考。
7. **日志与设置**：右上角切换中/英文界面并查看运行日志；左下角齿轮调整
   分析参数（如靶点过滤阈值）与 STRING 数据版本——官方发新版后在设置里
   改版本号并重建背景网络即可，无需更新程序。

所有中间结果自动缓存（见 `WY/config.py`）：四个靶点数据库的拉取缓存统一
60 天（上游发版周期约一季到两年半不等，不漏新版又把配额消耗压到最低），
整合表、PPI、富集等分析缓存默认 7 天。需要立即刷新时在管理页点「更新」
强制重新联网拉取。

## 命令行批量分析（CLI）

`WY/cli.py` 是 core 之上的薄封装，与网页共享全部缓存与设置：网页跑过的
分析 CLI 秒回，网页改过的阈值 CLI 自动沿用。CLI 面向批量出表，不绘图、
不对接，论文用图走网页端。

```bash
cd WY
source .venv/bin/activate   # Windows: .venv\Scripts\activate（启动脚本已建好环境）

# 单步命令：每个环节产出 CSV/JSON 到 --out 目录（默认 ./results，
# 可在 API.env 里用 CLI_RESULTS_DIR 配置固定输出目录）
python cli.py targets   --drug aspirin --disease "colorectal cancer"   # 四数据源靶点明细
python cli.py ppi       --genes genes.txt                              # STRING 网络 + MCC 核心靶点
python cli.py enrich    --genes common_targets.csv [--pvalue 0.05] [--sources GO:BP,KEGG]
python cli.py proximity --drug aspirin --disease "colorectal cancer"   # 邻近度 Z-score
python cli.py rwr       --seeds core_targets.csv                       # RWR 全网络排序

# 单对全流程：一个结果子目录 + summary.json（参数快照与各步状态，可追溯）
python cli.py analyze --drug aspirin --disease "colorectal cancer"

# 批量：pairs.csv 须含 drug、disease 两列，一行一对；每对一个子目录 + batch_summary.csv
python cli.py batch --pairs pairs.csv --out results
python cli.py batch --pairs pairs.csv --out results --resume   # 断点续跑
```

- **退出码**：全部成功 0，任一失败 1，失败原因直接打印（结果不会静默残缺）；
  单对失败不中断整批，`--resume` 跳过已成功的对、重试失败的对。
- **基因清单**三种写法对 `ppi`/`enrich` 的 `--genes` 与 `rwr` 的 `--seeds`
  通用：txt 一行一个基因（支持 `#` 注释）、csv（自动取 `gene`/`approvedSymbol`
  列）、逗号分隔列表。`targets`/`proximity`/`analyze`/`batch` 只收药名/病名，
  两组命令靠产出文件衔接（如 `enrich --genes common_targets.csv`）。
- proximity / RWR 依赖背景互作组：未构建时跳过并提示，不会自动触发大文件下载。
- 仅支持英文疾病名；批量建议顺序跑，不要同时开多个 batch 进程。
- `python cli.py serve` 可直接启动网页版服务器（自动打开浏览器）。

**装成全局命令 `npap`（可选）**：

```bash
cd WY
pip install -e ..        # 必须带 -e（可编辑安装）：缓存/日志/API.env 均相对
                         # WY/ 源码目录解析，普通安装会把它们写进 site-packages
npap --help              # 入口脚本装进当前环境的 Scripts\（Windows）或 bin/
```

装完后 `npap` 在已激活的环境里任意目录可用（三平台一致）；要免激活使用，
把该环境的 `Scripts\`（Windows）或 `bin/`（Linux/macOS）目录加入 PATH。
全局使用建议在 `API.env` 里配一次 `CLI_RESULTS_DIR`，否则不带 `--out` 时
结果会落到当前目录的 `results`。

## 配置说明

| 文件           | 用途                                                                                      |
|----------------|-------------------------------------------------------------------------------------------|
| `WY/API.env`   | 密钥与本机路径（`DisGeNET_API_KEY`、`VINA_DIR`、`OBABEL_DIR`、`CLI_RESULTS_DIR`），不入库 |
| `WY/config.py` | 数据源地址、限流/重试、熔断冷却阶梯、缓存 TTL、靶点过滤阈值、STRING 版本与下载并发        |

改行为看 `WY/config.py`，填密钥看 `WY/API.env`。

## 常见问题

- **DisGeNET 返回 429？** 程序已在学术配额（50 次/分）内本地限流，触发时稍等即可；
  配额为每天 25,000 次，正常单次分析仅消耗 1 次请求。本项目做了10%的余量配置，所以实际为45 次/分，每天20,000次。
  其他等级配额请自行调整 `WY/config.py` 的 `DISGENET_RATE_LIMIT_PER_MINUTE` / `_PER_DAY` / `_PER_MONTH`。
- **富集分析没结果？** g:Profiler 连续失败会自动熔断并改用 Enrichr；两源都不可用
  时页面提示"约 X 分钟后自动重试"。冷却阶梯见 `GPROFILER_CIRCUIT_COOLDOWN_LADDER`
  / `ENRICHR_CIRCUIT_COOLDOWN_LADDER`。
- **改了过滤阈值后结果没变？** 需同步递增 `TARGET_FILTER_VERSION`，旧整合缓存
  会自动失效重建。
- **背景互作组下载慢？** 国内镜像通常只同步最新版本，其它版本自动回退官方站；
  可调大 `WY/config.py` 的 `INTERACTOME_DOWNLOAD_CONCURRENCY`（默认 6）多段并行
  提速（官方站按连接限速），或先用镜像已同步的版本。
- **分子对接找不到 Vina / obabel？** 确认可执行文件已加入 PATH，或在 `WY/API.env`
  中将 `VINA_DIR` / `OBABEL_DIR` 指向相应目录（也可填可执行文件全路径）。
- **命令行结果没落在预期目录？** 落点优先级为 `--out` > `CLI_RESULTS_DIR`
  （`WY/API.env` 或同名系统环境变量）> 当前目录的 `results`；`--out` 可单次覆盖。

## 项目结构

```
WY/
├── api.py                  # FastAPI 入口：/api 端点 + 托管 web/dist 构建产物
├── cli.py                  # CLI 入口：serve / targets / ppi / enrich / proximity / rwr / analyze / batch
├── config.py               # 全部可调参数（路径、限流/重试、缓存 TTL、过滤阈值、STRING 版本）
├── env_discovery.py        # Vina / Open Babel 可执行文件探测
├── API.env                 # 密钥与本机路径（不入库；设置页首次保存后自动生成）
├── core/                   # 领域与编排层，不依赖任何 UI 框架
│   ├── targets / sources / chembl / dgidb / disgenet / opentargets
│   │                       #   四数据源抓取与靶点整合降噪
│   ├── ppi / enrichment    #   STRING 网络、GO/KEGG 富集（双源熔断）
│   ├── interactome / network_algorithms    #   背景互作组、邻近度 / RWR
│   ├── docking / docking_session           #   AutoDock Vina 对接
│   ├── geo / expression / ml_screening / ml_pipeline / ml_results
│   │                       #   GEO 表达谱 + ML 共识筛选
│   ├── management / jobs / analysis / files  #  数据管理编排、后台任务、结果文件
│   └── http_client / caching / logging_setup / text_utils  # 基础设施
├── web/                    # TS + Vue 3 前端源码（Vite 构建，产物 dist/ 由 api.py 托管）
├── i18n/                   # 中英文案词典与取词 API
├── tools/                  # 独立工具：build_interactome / geo_ml_screen / net_block / export_i18n_json
├── benchmark/              # 基准验证：LIT-PCBA / 跨数据集稳健性 / ML 规模 / 恢复性
└── Static/                 # 运行时产生的缓存与分析结果
```

依赖严格单向：`web(TS+Vue) / api.py → core/ → i18n、config`。`core/` 不依赖
任何 UI 框架，可在无界面环境下导入、复用与测试；`cli.py` 只依赖 `core/` 与
`config`，与网页共享缓存与设置。

## 开发与测试

```bash
pytest                            # 全套测试，纯本地，外部请求全部打桩
pytest --cov --cov-report=html    # 覆盖率报告（生成 htmlcov/，浏览器打开逐行红绿图）
pytest -p WY.tools.net_block      # 封网审计：拦截非回环连接，证明零外网依赖
```

- 覆盖率口径见 `pyproject.toml`：统计 `WY` 包（不含 `tools/`、
  `benchmark/`、`web/`、`Static/`），低于 75% 直接判失败，本地与 CI 同一口径。
- ML 相关测试默认瘦身档，设 `NPAP_ML_TEST_FULL=1` 切全参。
- CI（`.github/workflows/ci.yml`）：每次 push / PR 自动从零装依赖并跑全套测试；
  每周一另跑一次全参定时档，兜住依赖升级带来的统计行为漂移。
- 外部数据源（DGIdb / DisGeNET / ChEMBL / Open Targets / STRING / g:Profiler /
  Enrichr / GEO）一律打桩，写盘由 `tests/conftest.py` 沙箱接管：不消耗 API
  配额、不触碰真实 `Static/` 缓存。
- 背景互作组构建：优先在数据管理页点「一键构建」，或运行
  `WY/tools/build_interactome.py`。
- 分子对接基准验证：`WY/benchmark/lit_pcba_benchmark.py`（其余基准脚本同目录）。

## 免责声明

本项目仅供科研用途，分析结果请以官方数据库与实验验证为准。
使用数据时请遵循各数据源（STRING、DisGeNET、ChEMBL、DGIdb、Open Targets、GEO）
的引用要求。第三方组件许可清单见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
若本项目有任何问题，请联系项目作者。
`NPAP` 为本项目自用缩写，与任何第三方机构、项目或商标无关。