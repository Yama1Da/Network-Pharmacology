简体中文 | [EN](README_EN.md)

# 网络药理学分析平台

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

以下要求针对 conda / pip 两种本地方式；用 Docker 部署时，宿主机只需 Docker 与 DisGeNET 密钥（见下方「Docker 部署」）。

- Python >= 3.13, < 3.14（开发基准 3.13.x；依赖钉版尚未提供 3.14 的预编译轮子）
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) 可执行文件（按平台安装方式见「快速开始」第 3 步；在 `WY/API.env`（请将API.env.example重命名为API.env）中指定 `VINA_DIR`，或将其所在目录加入 PATH）
- [Open Babel](https://openbabel.org/) 的 `obabel` 命令行（分子对接受体转换以子进程调用；按平台安装方式见「快速开始」第 3 步；可指定 `OBABEL_DIR`，或将其所在目录加入 PATH）
- [DisGeNET](https://www.disgenet.org/) 学术账号 API Key（免费申请）/ 其他级别账号 API Key

## 快速开始

```bash
# 1. 安装依赖（二选一）

# 方式 A（推荐）：conda 一条命令建依赖环境（obabel 由此安装，三平台通用；
#               密钥与本机路径不属于环境，仍需下面第 2、3 步配置）
cd WY
conda env create -f environment.yml
conda activate npy

# 方式 B：pip + venv
cd WY
python -m venv .venv
.venv\Scripts\activate        # Windows（Linux/macOS: source .venv/bin/activate）
pip install -r requirements.txt
# 受体转换还需要 obabel 命令行（子进程调用，非 Python 绑定）：Windows 去官方网站装 Open Babel 官方安装包、
# Linux 用 apt/dnf 安装 openbabel，或 conda install -c conda-forge openbabel；
# 保证 obabel 在 PATH，或在 API.env 里用 OBABEL_DIR 指向其所在目录

# 2. 配置 API.env（conda / pip 两种方式完全相同：密钥是运行时配置，不属于任何环境）
复制 API.env.example 为 API.env（保留模板便于日后接收更新），编辑并填入 DisGeNET_API_KEY

# 3. 安装 AutoDock Vina 与 Open Babel（vina 在 conda-forge 无 win-64 构建；obabel 供受体转换）
Windows：到 https://github.com/ccsb-scripps/AutoDock-Vina/releases 下载
#          vina_1.2.7_win.exe（可连同 vina_split_1.2.7_win.exe 一起）放入任意目录，
#          在 API.env 中将 VINA_DIR 指向该目录（或将该目录加入 PATH）；
#          不装 vina_split 亦可，程序已内置 python 分割兜底。
#          Open Babel 到 https://openbabel.org/ 下载 Windows 安装包（含 obabel.exe），
#          勾选加入 PATH 后即可；否则在 API.env 中将 OBABEL_DIR 指向其所在目录。
Linux/macOS：conda install -c conda-forge vina   # 同时安装 vina 与 vina_split
#          obabel：conda install -c conda-forge openbabel（或 apt/dnf 安装 openbabel）

# 4. 启动
streamlit run App.py

# 5. （可选）把命令行批量分析装成全局命令 npap（详见「命令行批量分析（CLI）」一节）
pip install -e ..        # 在 WY 目录下执行；必须带 -e（可编辑安装）
npap --help              # 在已激活的环境中任意目录可用
```

## Docker 部署（Linux 服务器 / 对外分发）

面向 Linux 服务器部署、对外分发与答辩演示；Windows 科研用户建议用上面的 conda 方式。
宿主机只需 Docker（含 compose 插件）与 DisGeNET 密钥，无需装 Python、Vina、Open Babel——
这些全部在镜像内（`vina` / `obabel` 已在 PATH，容器内无需 `VINA_DIR` / `OBABEL_DIR`）。

```bash
# 1. 准备钥匙：密钥不在镜像里，运行时由 compose 注入
cp WY/API.env.example API.env    # 编辑填入 DisGeNET_API_KEY（空 key 也能启动，仅 DisGeNET 靶点表为空）

# 2. 准备缓存目录（容器以 uid 1000 运行；须先建好——compose 自动创建的目录属主是 root）
mkdir -p Static && sudo chown 1000:1000 Static

# 3. 构建并启动（首次构建较慢，镜像约 1.2 GB；之后改代码重建只需秒级到分钟级）
docker compose up -d --build
# 浏览器打开 http://<服务器IP>:8501；日志: docker logs -f network-pharmacy
```

- **密钥模型**：镜像 = 程序，`API.env` = 钥匙。`API.env` 与 `Static` 数据被 `.dockerignore`
  挡在构建上下文外，不会进镜像；compose 的 `env_file` 在运行时逐行注入容器环境变量。
- **缓存持久化**：`./Static` 挂载进容器，重建/升级容器缓存不丢。
- **对外分发**：默认 `docker compose up -d --build` 在本机构建；也可把镜像推到 registry
  （Docker Hub / 阿里云 ACR 等），把 compose.yaml 里的 `build: .` 注释掉、`image:` 换成
  registry 路径，对方免构建直接 `docker compose up -d`。弱网环境可用
  `docker save | gzip` 导出 tar 包离线搬运。
- **容器内跑 CLI**（镜像内未装 `npap` 全局命令，统一用 `python cli.py` 调用；结果经
  `--out` 落到挂载目录，宿主机直接可取；批量时把 `pairs.csv` 放进宿主机 `./Static/`，
  容器内以 `/app/WY/Static/pairs.csv` 引用）：

```bash
docker compose run --rm app python cli.py analyze --drug aspirin --disease "colorectal cancer" \
    --out /app/WY/Static/cli-results
```

- **升级版本**：`git pull && docker compose up -d --build`（依赖未变时只重建代码层，很快）。

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

## 命令行批量分析（CLI）

`WY/cli.py` 是 core 之上的薄封装：与网页共享全部缓存与设置（网页跑过的对 CLI 秒回；
网页齿轮改过的阈值 CLI 自动沿用），面向「湿实验 + 文献对照」的批量出表。
不绘图、不对接——论文用图走网页端（图均可下载），CLI 保持纯数据输出。

```bash
cd WY
conda activate npy        # 或激活.venv

# 单步命令：每个环节产出 CSV/JSON 到 --out 目录（默认 ./results；
# 可在 API.env 里用 CLI_RESULTS_DIR 配置一次永久输出目录）
python cli.py targets   --drug aspirin --disease "colorectal cancer"   # 四数据源靶点明细
python cli.py ppi       --genes genes.txt                              # STRING 网络 + MCC 核心靶点
python cli.py enrich    --genes common_targets.csv [--pvalue 0.05] [--sources GO:BP,KEGG]
python cli.py proximity --drug aspirin --disease "colorectal cancer"   # 邻近度 Z-score
python cli.py rwr       --seeds core_targets.csv                       # RWR 全网络排序

# 单对全流程：一个结果子目录 + summary.json（参数快照与各步状态，结果可追溯）
python cli.py analyze --drug aspirin --disease "colorectal cancer"

# 批量跑批：pairs.csv 表头须含 drug、disease 两列（大小写不敏感，Excel 里的
# No. 序号列等多余列自动忽略；Excel 编辑后另存即可），
# 行 = 一对，一对一配对；每对一个子目录 + batch_summary.csv 汇总
python cli.py batch --pairs D:\cases\pairs.csv --out D:\cases\results
python cli.py batch --pairs D:\cases\pairs.csv --out D:\cases\results --resume   # 断点续跑
```

- **退出码**：全部成功 0，任一失败 1（CI/自检脚本可直接用）；失败会在终端给出具体
  原因（富集双源不可用、疾病名无返回、背景互作组未构建等），结果不会静默残缺。
- **单对失败不中断整批**；`--resume` 跳过已成功的对、自动重试失败的对。
- 基因清单三种写法对 `ppi`/`enrich` 的 `--genes` 与 `rwr` 的 `--seeds` **通用**，任选其一：
  txt 一行一个基因（支持 `#` 注释）、csv（自动取 `gene`/`approvedSymbol` 列）、逗号分隔
  字面列表。
  名称类命令（`targets`/`proximity`/`analyze`/`batch`）只收药名/病名。
  两组命令靠产出文件衔接（如 `enrich --genes <对目录>\common_targets.csv`）。
- proximity / RWR 依赖背景互作组：未构建时跳过并提示，不会自动触发大文件下载
  （构建走网页「一键构建」或 `scripts/build_interactome.py`）。
- 与网页一样仅支持英文疾病名；批量时建议顺序跑，不要同时开多个 batch 进程。

**装成全局命令 `npap`（可选；开发期 `python cli.py` 即可用；Windows/Linux/macOS 通用）**：

```bash
cd WY
pip install -e ..        # conda 用户先 conda activate npy；必须带 -e（可编辑安装），见下
npap --help              # 入口脚本装进当前环境的 Scripts\（Windows）或 bin/（Linux/macOS）
```

- **必须用可编辑安装（`-e`）**：项目按"就地运行"设计，缓存/日志/`API.env` 均相对
  `WY/` 源码目录解析；不带 `-e` 的普通安装会把它们写进 site-packages。
- 装完后 `npap` 在**已激活的环境**里任意目录可用（三平台一致）；要在不激活环境的终端
  直接用，需把该环境的 `Scripts\`（Windows）或 `bin/`（Linux/macOS）目录加入 PATH。
- Windows 专属的免安装过渡：复制 `npap.bat.example` 为 `npap.bat`，把所在目录加入 PATH
  即可（含本机路径，与 `API.env` 同类不入库；标准布局零修改可用，conda 等其他布局把其中
  的 python 路径改成绝对路径）。Linux/macOS 免装环境可用 Docker 直接跑（见下方「Docker 部署」）。
- 全局使用建议先在 `API.env` 里配一次 `CLI_RESULTS_DIR`（结果输出目录）：否则不带 `--out`
  时结果会落到**当前目录**的 `results` 目录。

## 配置说明

| 文件 | 用途 |
|---|---|
| `WY/API.env` | 密钥与本机路径（`DisGeNET_API_KEY`、`VINA_DIR`、`OBABEL_DIR`、`CLI_RESULTS_DIR`），不入库 |
| `WY/config.py` | 数据源地址、限流/重试、熔断冷却阶梯、缓存 TTL、靶点过滤阈值、STRING 版本与下载并发 |

> 改行为看 `WY/config.py`，填密钥看 `WY/API.env`：`WY/API.env` 只放不提交到 git 的敏感信息
> 与本机路径，其余请求控制、缓存有效期、靶点降噪过滤阈值集中在 `WY/config.py`。

## 项目结构

```
WY/
├── App.py                  # 首页：发起分析、靶点概览、全局设置入口
├── cli.py                  # 命令行批量分析入口（targets/ppi/enrich/proximity/rwr/analyze/batch）
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
不依赖任何 UI 框架（Streamlit），可在无界面环境下导入、复用与测试；
`cli.py` 同样只依赖 `core/` 与 `config`，与网页共享缓存与设置。

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
- **分子对接找不到 Vina / obabel？** 确认可执行文件已加入 PATH，或在 `WY/API.env` 中将
  `VINA_DIR` / `OBABEL_DIR` 指向相应目录（也可填可执行文件全路径）。受体转换用的 `obabel`
  来自 Open Babel，安装方式见「快速开始」第 3 步。
- **命令行结果没落在预期目录？** 不带 `--out` 时的落点优先级为：`--out` 参数 >
  `CLI_RESULTS_DIR`（`WY/API.env` 或同名系统环境变量）> 当前目录的 `results`；
  `--out` 仍可单次覆盖，无需删除已配置的 `CLI_RESULTS_DIR`。
- **Docker 起不来或页面打不开？** 先看 `docker compose ps` 与 `docker logs network-pharmacy`
  定位；远程服务器上浏览器应访问 `http://<服务器IP>:8501` 而非 localhost，
  并确认防火墙/安全组已放行 8501。

## 免责声明

本项目仅供科研用途，分析结果请以官方数据库与实验验证为准。
使用数据时请遵循各数据源（STRING、DisGeNET、ChEMBL、DGIdb、Open Targets）的引用要求。