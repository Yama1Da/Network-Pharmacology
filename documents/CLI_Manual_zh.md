# NPAP 命令行(CLI)操作手册
[简体中文](#中文版) | [English](CLI_Manual_en.md)

---

## 中文版

### 0. 前置条件

1. 已用启动脚本(`start_npap.bat` / `start_npap.sh`)装好环境,虚拟环境在 `WY/.venv`:
   ```bash
   cd WY
   source .venv/bin/activate        # Windows: .venv\Scripts\activate
   python cli.py --help             # 或装成全局命令后直接 npap --help
   ```
2. **背景互作组**(proximity / rwr / modules / synergy 的前置):在网页
   「数据管理」页点「一键构建」,或运行 `WY/tools/build_interactome.py`。
   未构建时这些命令会警告并跳过,不会自动触发大文件下载。
3. **DisGeNET API Key**(可选):不填则疾病靶点表少一个来源,其余不受影响。
4. 疾病名**仅支持英文**(如 `Type 2 diabetes mellitus`)。

### 1. 子命令总览

| 子命令      | 用途                                                    | 典型耗时*                   |
|-------------|---------------------------------------------------------|-----------------------------|
| `serve`     | 启动网页版服务器(自动开浏览器)                          | 秒级                        |
| `targets`   | 四数据源(ChEMBL/DGIdb/Open Targets/DisGeNET)靶点明细    | ~1 分钟(有缓存则秒回)       |
| `ppi`       | STRING 网络 + 拓扑指标 + MCC 核心靶点                   | 秒~分钟                     |
| `enrich`    | GO/KEGG 富集(g:Profiler,失败自动切 Enrichr)             | ~30 秒                      |
| `proximity` | 药物-疾病靶点模块邻近度 Z-score                         | ~2 分钟(200 置换)           |
| `rwr`       | 种子基因 RWR 全网络排序                                 | ~1 分钟                     |
| `analyze`   | 单对全流程(靶点→交集→PPI→富集→邻近度→RWR)               | 3~5 分钟                    |
| `batch`     | 多对全流程 + 断点续跑 + 批量汇总表                      | 对数 × analyze              |
| `modules`   | 疾病模块三方法(RWR/DIAMOnD/近邻)+ Borda 共识 + 置换检验 | 首跑分钟级(先小 perms 试跑) |
| `synergy`   | 网络分离度协同预测(单对或批量)                          | 分钟级/对                   |

\* 首跑含联网取数;缓存命中(默认 TTL 7 天,靶点源 60 天)后大幅缩短。

### 2. 通用约定

- **输出目录优先级**:`--out` > `CLI_RESULTS_DIR`(配置在 `WY/API.env` 或同名
  环境变量)> 当前目录 `./results`。
- **退出码**:全部成功 0;任一失败 1(原因直接打印,结果不静默残缺)。
- **基因清单三种写法**(`--genes`/`--seeds` 通用):
  ① `.txt` 一行一个基因(支持 `#` 注释);② `.csv`(自动取 `gene`/
  `approvedSymbol` 列);③ 逗号分隔字符串,如 `"TP53,EGFR"`。
  `targets`/`proximity`/`analyze`/`batch` 只收药名/病名;两组命令靠产出文件衔接。
- **与网页同源**:同一对参数网页/CLI 共享缓存与设置页阈值。
- **批量纪律**:顺序跑,不要同时开多个 batch 进程(缓存写竞争)。

### 3. 逐命令说明

#### 3.1 `targets` — 四库靶点明细

```bash
python cli.py targets --drug aspirin --disease "colorectal cancer"
python cli.py targets --drug aspirin                # 只取药物侧
```
产物:`<drug>_targets.csv` / `<disease>_targets.csv` / `<drug>_<disease>_common_targets.csv`。
注意:DisGeNET 报空先检查 `WY/API.env` 里的 key;数据源失败按空表处理并警告,
**退出码 1**(明细不完整,文献对照不可用)。

#### 3.2 `ppi` — STRING 网络与核心靶点

```bash
python cli.py ppi --genes common_targets.csv
```
产物:边表 / 指标表(度数、中介中心性等)/ MCC 核心靶点表。基因数超过
STRING 上限自动截断并警告;与网页"手动基因模式"共享边表缓存。

#### 3.3 `enrich` — GO/KEGG 富集

```bash
python cli.py enrich --genes common_targets.csv --pvalue 0.05 --sources GO:BP,KEGG
```
产物:富集长表 CSV。g:Profiler 连续失败自动熔断切换 Enrichr,双源均不可用
时提示"约 X 分钟后自动重试"。

#### 3.4 `proximity` — 邻近度 Z-score

```bash
python cli.py proximity --drug aspirin --disease "colorectal cancer"
```
产物:`proximity.json`(Z-score、经验 p 值、度匹配零模型参数)。未建互作组
时警告跳过(退出码 1)。

#### 3.5 `rwr` — RWR 全网络排序

```bash
python cli.py rwr --seeds core_targets.csv
python cli.py rwr --seeds "TP53,EGFR,MTOR"
```
产物:`rwr.csv`(gene, score, is_seed;得分降序即优先级)。

#### 3.6 `analyze` — 单对全流程

```bash
python cli.py analyze --drug aspirin --disease "colorectal cancer"
```
产物:结果子目录 + `summary.json`(**参数快照 + 每步状态/耗时**,可追溯)。
单步失败不阻断后续;结束打印未成功的步骤名单。

#### 3.7 `batch` — 批量全流程

```bash
python cli.py batch --pairs pairs.csv --out results
python cli.py batch --pairs pairs.csv --out results --resume   # 断点续跑
```
`pairs.csv` 须含 `drug`、`disease` 两列(一行一对)。产物:每对一个子目录 +
`batch_summary.csv`(状态/耗时/失败原因,一行一对对得上账)。`--resume`
跳过已成功的对、重试失败的对;任一失败退出码 1。

#### 3.8 `modules` — 疾病模块三方法一致性

```bash
python cli.py modules --genes disease_targets.csv                     # 默认 N=100,跑置换检验
python cli.py modules --genes disease_targets.csv --top-n 20 --n-perms 10   # 小参数试跑
python cli.py modules --genes disease_targets.csv --n-perms 0        # 跳过显著性
```
参数:`--top-n` 每方法扩张规模 N(默认 100,三方法共用);`--k` 一致性阈值
(默认 2:被 ≥2 种方法捞入才算一致);`--n-perms` 置换轮数(默认 200;
**DIAMOnD 零分布真实图分钟级起,先小轮数试跑**)。
产物:`consensus.csv`(三方法秩 + Borda 共识分 + consistent 标记)/
`modules.csv`(三方法成员长表)/ `significance.csv`(各方法置换 p 值)/
`summary.json`(参数快照;LCC 外种子进 `seeds_out_of_lcc` 不静默丢)。

#### 3.9 `synergy` — 网络分离度协同预测

```bash
# 单对
python cli.py synergy --a Paeonol --b Luteolin --disease "breast carcinoma"
# 无疾病锚(纯分离度)也可
python cli.py synergy --a Paeonol --b Luteolin
# 批量:CSV 须含 drug_a 列;drug_b/disease 列可选(缺省用 --b/--disease 固定值)
python cli.py synergy --input pairs.csv --out results
python cli.py synergy --input pairs.csv --b melatonin --top-n 50 --n-perm 200
```
行为:每对产出契约 JSON 归档进共享 `synergy-results/` 库(与网页结果管理页
互见);批量另出 `synergy_batch_summary.csv`(**含 E6 可达性五列**:
r_tier_a / r_tier_b / r_high_conf / r_sort_key / r_sort_key_note——排序比较
直接在 Excel 里按 r_sort_key 降序)。`--top-n` 靶点切片(0=全量,默认 100);
`--n-perm` 置换轮数(≥10,默认 200)。逐对 CSV 明细在各自运行目录
(`reachability_<药>.csv`)。CLI 不绘图;p 值不稳定时结果行带 `power_limited`。

#### 3.10 `serve` — 启动网页版

```bash
python cli.py serve                    # 自动打开浏览器
python cli.py serve --port 8000 --no-browser
```

### 4. 端到端教程

**教程 A:单药物完整网络药理学分析**
```bash
python cli.py analyze --drug aspirin --disease "colorectal cancer"
# 想分步检查就把 analyze 拆成 targets → ppi/enrich → proximity/rwr
```

**教程 B:批量筛查 + 断点续跑**
```bash
# 1) 准备 pairs.csv(两列:drug,disease)
# 2) 跑批(中断后加 --resume 续)
python cli.py batch --pairs pairs.csv --out results --resume
# 3) 检查 batch_summary.csv 里 failed 行的原因
```

**教程 C:协同预测批量 + 可达性排序**
```bash
# 1) 准备 pairs.csv(列:drug_a, drug_b, disease;disease 可留空 = 无疾病锚)
# 2) 跑批
python cli.py synergy --input pairs.csv --out results
# 3) 打开 synergy_batch_summary.csv,按 r_sort_key 降序 = E6 全局排序
#    (r_high_conf=True 的行 = 两药都有暴露证据,可信度最高)
```

### 5. 常见错误对照

| 现象                                     | 原因与处理                                                                                     |
|------------------------------------------|------------------------------------------------------------------------------------------------|
| `background interactome not built`       | proximity/rwr/modules/synergy 前置未满足:数据管理页一键构建或跑 `tools/build_interactome.py`   |
| `cannot write 'xxx.csv': file is locked` | Excel 占用输出文件,关掉重跑                                                                    |
| DisGeNET 靶点为空 / 报 429               | 检查 `WY/API.env` 密钥;程序已限流(45 次/分),触发稍等                                           |
| `no genes parsed from --genes`           | 文件编码/列名问题:csv 需含 `gene` 或 `approvedSymbol` 列                                       |
| 结果没落在预期目录                       | 落点优先级 `--out` > `CLI_RESULTS_DIR` > `./results`;建议在 `API.env` 配一次 `CLI_RESULTS_DIR` |
| synergy 批量某对 `synergy:no_targets_*`  | 四库查无该药/病靶点(多为拼写);其余对不受影响,退出码 1                                          |
| 改了阈值结果没变                         | 需同步递增 `WY/config.py` 对应 `*_VERSION`,旧缓存自动失效重建                                  |

### 6. 全局命令 `npap`(可选)

```bash
cd WY && pip install -e ..    # 必须带 -e:缓存/日志相对 WY/ 源码目录解析
npap --help                   # 之后任意目录可用(激活环境内)
```
建议在 `API.env` 配一次 `CLI_RESULTS_DIR`,否则不带 `--out` 时结果落到
当前目录的 `results`。
