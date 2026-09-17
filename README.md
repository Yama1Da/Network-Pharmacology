简体中文 | [EN](README_EN.md)

# 网络药理学分析平台（Network Pharmacology Platform）

基于 Streamlit 的一站式网络药理学分析工具：输入药物与疾病名称，自动完成
「靶点获取 → 靶点整合 → PPI 网络构建 → 富集分析 → 分子对接」的完整流程。


## 功能特性

- **疾病/药物靶点获取**：接入 DisGeNET、DGIdb、Open Targets、ChEMBL 四个数据源，
  自动整合去重，并过滤非特异性结合蛋白（白蛋白、CYP 代谢酶等）降噪，降噪条件可以自己调参。
- **PPI 网络**：基于 STRING 构建蛋白互作网络（默认 v12.0，设置面板可切换数据版本），支持多种可视化风格展示
- **背景互作组**：STRING 全网络本地缓存，支撑邻近度 Z-score / RWR 靶点优先级分析；
  数据管理页一键构建，镜像/官方双源自动回退 + 分段并发下载
- **富集分析**：GO / KEGG 等通路富集，g:Profiler 失败时自动切换 Enrichr；
  两源均带熔断降级（冷却阶梯逐次加长、到期自动探测恢复），双源都不可用时页面提示
  "约 X 分钟后自动重试"
- **分子对接**：集成 AutoDock Vina，RDKit + Meeko 自动完成配体准备，内置 3D 结果查看 + 聚焦口袋（虚化周边）
- **工程化**：请求级缓存、DisGeNET 本地限流（学术配额内）、富集数据源熔断降级与自动恢复、
  中英双语界面、每日滚动日志、靶点/疾病数据管理页面

## 环境要求

- Python >= 3.13, < 3.14（开发基准 3.13.x；依赖钉版尚未提供 3.14 的预编译轮子）
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) 可执行文件（按平台安装方式见「快速开始」第 3 步；在 `WY/API.env`（请将API.env.example重命名为API.env）中指定 `VINA_DIR`，或将其所在目录加入 PATH）
- [DisGeNET](https://www.disgenet.org/) 学术账号 API Key（免费申请）/ 其他级别账号 API Key

## 快速开始

```bash
# 1. 安装依赖（二选一）

# 方式 A（推荐）：conda 一条命令建依赖环境（openbabel 由此安装，三平台通用；
#               密钥与本机路径不属于环境，仍需下面第 2、3 步配置）
cd WY
conda env create -f environment.yml
conda activate npy

# 方式 B：pip + venv
cd WY
python -m venv .venv
.venv\Scripts\activate        # Windows（Linux/macOS: source .venv/bin/activate）
pip install -r requirements.txt
# openbabel 不在 requirements.txt（pip 安装在缺编译环境的平台常失败）：
# 应用可正常启动，仅分子对接的受体转换需要它（懒加载），届时补装：
# conda install -c conda-forge openbabel    （或 pip install openbabel）

# 2. 配置 API.env（conda / pip 两种方式完全相同：密钥是运行时配置，不属于任何环境）
复制 API.env.example 为 API.env（保留模板便于日后接收更新），编辑并填入 DisGeNET_API_KEY

# 3. 安装 AutoDock Vina（conda-forge 的 vina 无 win-64 构建，按平台二选一）
# Windows：到 https://github.com/ccsb-scripps/AutoDock-Vina/releases 下载
#          vina_1.2.7_windows_x86_64.exe 放入任意目录，在 API.env 中将 VINA_DIR
#          指向该目录（或将该目录加入 PATH）；vina_split 无官方 Windows 版，
#          程序已内置 python 分割兜底，可不装。
# Linux/macOS：conda install -c conda-forge vina   # 同时安装 vina 与 vina_split

# 4. 启动
streamlit run App.py
```

## 使用流程

1. **首页发起分析**：输入药物名与疾病名（默认示例为 metformin 与 Type 2 diabetes mellitus）,现暂时不支持中文疾病名输入，
   点击「开始分析」。程序会并发请求 ChEMBL / DGIdb / Open Targets / DisGeNET 四个数据源。
2. **查看靶点概览**：分析完成后首页展示药物靶点数、疾病靶点数、交集靶点数与韦恩图，
   展开折叠面板可查看各数据源的靶点明细表（按来源数与打分排序）。
3. **数据管理**：通过首页左上角的「数据管理」菜单进入各信息页面，
   查看、增删或修正抓取到的疾病靶点与药物靶点；同一页面可一键构建背景互作组
   （邻近度 / RWR 分析的前置）。
4. **PPI 网络与富集**：从菜单进入「数据网络」页面，基于整合后的靶点构建 STRING
   互作网络（可选多种可视化风格），并执行 GO / KEGG 富集分析。
5. **分子对接**：从菜单进入「分子对接」页面，选择核心靶点进行 AutoDock Vina 对接，
   在线查看 3D 结合构象与聚焦口袋效果。
6. **日志与设置**：右上角可切换中/英文界面并查看运行日志；左下角齿轮可调整全局设置
   （如靶点过滤阈值）与 STRING 数据版本（官方发布新版本后，在设置里输入版本号并重建
   背景网络即可，无需更新程序），设置变更后整合缓存自动失效重建。

所有中间结果自动缓存（默认 7 天，见 `WY/config.py`），重复分析秒级返回。（由于DGIdb 不定期更新，所以设置缓存 7 天以上。）

## 配置说明

| 文件 | 用途 |
|---|---|
| `WY/API.env` | 密钥与可执行文件路径（`DisGeNET_API_KEY`、`VINA_DIR`），不入库 |
| `WY/config.py` | 数据源地址、限流/重试、熔断冷却阶梯、缓存 TTL、靶点过滤阈值、STRING 版本与下载并发 |

> 改行为看 `WY/config.py`，填密钥看 `WY/API.env`：`WY/API.env` 只放不提交到 git 的敏感信息
> 与本机路径，其余请求控制、缓存有效期、靶点降噪过滤阈值集中在 `WY/config.py`。

## 项目结构

```
WY/
├── App.py                  # 首页：发起分析、靶点概览、全局设置入口
├── config.py               # 全部可调参数（路径、限流/重试、缓存 TTL、过滤阈值）
├── API.env                 # 密钥与本机路径（不入库，模板见 API.env.example）
├── pages/                  # Streamlit 子页面（靶点管理、数据网络、富集、对接、日志等）
├── core/                   # 领域与编排层：数据源抓取、PPI、富集、对接、缓存文件管理
│   ├── management.py       #   靶点数据文件刷新编排、设置持久化、缓存文件名可读化
│   ├── targets.py          #   多数据源靶点获取与整合
│   ├── ppi.py / enrichment.py / docking.py / interactome.py / network_algorithms.py
│   └── http_client.py / caching.py / text_utils.py / logging_setup.py  # 基础设施
├── ui/                     # 界面组件层
│   ├── widgets.py          #   跨页面复用小组件（文件表、批量操作栏、搜索框等）
│   └── source_management.py  # 药物/疾病靶点管理页整页渲染
├── i18n/                   # 中英文案词典与取词 API
├── scripts/                # 离线自检脚本（见「开发与自检」）
└── Static/                 # 运行时产生的缓存与分析结果
```

依赖严格单向：`pages / App.py → ui/ → core/ → i18n、config`。`core/` 与 `i18n/`
不依赖任何 UI 框架（Streamlit），可在无界面环境下导入、复用与测试。

## 开发与自检

```bash
cd WY
python scripts/selfcheck_docking.py   # 分子对接模块自检
python scripts/selfcheck_network.py   # 网络算法与富集模块自检
python scripts/build_interactome.py   # 如果不知道背景互作组是否已下载，可运行此脚本下载
```

前两个自检脚本均为纯本地运行（外部请求全部打桩），不消耗数据源配额；
修改 `core/` 后建议跑一遍作为回归。
**build_interactome.py**脚本建议在数据管理页面点击「一键构建背景互作组」运行，以确保背景互作组已下载。


## 常见问题

- **DisGeNET 返回 429？** 程序已在学术配额（50 次/分）内做了本地限流，仍触发时稍等即可；
  配额为每天 25,000 次，正常单次分析仅消耗 1 次请求。如有其他等级配额，请自行更改 `WY/config.py`中的
  `DISGENET_RATE_LIMIT_PER_MINUTE` `DISGENET_RATE_LIMIT_PER_DAY` `DISGENET_RATE_LIMIT_PER_MONTH`。
- **富集分析没结果？** g:Profiler 连续请求失败会自动熔断并改用 Enrichr；两源都不可用时页面
  会提示"约 X 分钟后自动重试"。冷却阶梯（30 秒起、逐次加长至 30 分钟上限）见
  `GPROFILER_CIRCUIT_COOLDOWN_LADDER` / `ENRICHR_CIRCUIT_COOLDOWN_LADDER`。
- **改了过滤阈值后结果没变？** 靶点过滤阈值变更需同步递增 `TARGET_FILTER_VERSION`，
  旧整合缓存会自动失效重建。
- **背景互作组下载慢？** 国内镜像通常只同步最新版本（实测目前仅 v12.0），其它版本会自动
  回退官方下载站（单连接较慢）。可将 `WY/config.py` 的 `INTERACTOME_DOWNLOAD_CONCURRENCY`
  （默认 6）调大提速——官方站按连接限速，多段并行近似线性提速；或先用镜像已同步的版本。
- **分子对接找不到 Vina？** 确认 `vina` 已加入 PATH 或环境变量存在 `VINA_DIR` ，或在 `WY/API.env` 中将 `VINA_DIR`
  指向 Vina 所在目录（也可填可执行文件全路径）。

## 免责声明

本项目仅供科研用途，分析结果请以官方数据库与实验验证为准。
使用数据时请遵循各数据源（STRING、DisGeNET、ChEMBL、DGIdb、Open Targets）的引用要求。