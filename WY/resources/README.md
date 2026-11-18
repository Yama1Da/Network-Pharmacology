# WY/resources —— 人工精选数据（源数据，进仓库）

本目录放**人工整理、需要版本化追溯**的源数据，与 `WY/Static/`（运行时缓存、不进仓库）
严格区分：这里的每一行都要能说清"谁、何时、依据哪份说明书"。

## dailymed_primary_targets.csv

二期 1.10 的说明书机制靶点真值表，供 1.7 可达性恢复基准使用：
**真值 = FDA 说明书明确标注的机制靶点基因**；说明书写"β 受体激动剂"这类类靶点
而不指明亚型的，填"待定"，不计入真值（宁缺毋滥）。

### 列定义

| 列                        | 含义                                                                               |
|---------------------------|------------------------------------------------------------------------------------|
| `drug`                    | 药物通用名（与平台检索名一致，可追溯）                                             |
| `label_mechanism_targets` | **人工确认**的机制靶点基因符号，分号分隔；无法判定写 `待定`                        |
| `ot_moa_prefill`          | Open Targets MOA 预填候选（自动生成，**不是真值**）                                |
| `spl_pharm_class`         | FDA 说明书药理类别原文（EPC / MoA / CS，分号分隔）                                 |
| `label_section`           | 来源定位：`<SPL setid>#43679-0`（43679-0 即 12.1 Mechanism of Action 的 LOINC 码） |
| `decision`                | `y`=确认、`t`=待定、`x`=否决、空=未审查                                            |
| `curator`                 | 填表人（改一条也要可追溯）                                                         |
| `date`                    | 审查日期 `YYYY-MM-DD`                                                              |
| `notes`                   | 机制原文摘录、否决原因等                                                           |

### 工具用法

```bash
# 1) 预填（三级：OT MOA → openFDA 药理类别 → DailyMed MoA 原文）
python WY/tools/dailymed_curate.py prefill --drugs drugs.txt --curator 你的名字

# 2) 逐条审查（y/t/x，每键立即写盘，可 Ctrl+C 后续跑）
python WY/tools/dailymed_curate.py review --curator 你的名字

# 3) 校验（格式硬校验；--online 才检查基因符号是否存在）
python WY/tools/dailymed_curate.py validate --online
```

离线语义：有缓存用缓存；无缓存时**警告并跳过**该药/该校验，不静默失败。

## cancermine_cancer_genes.csv

二期 1.9 泛癌对照基因集的 CancerMine 部分（CC0）。**不手改**，由
`WY/tools/build_pancancer_module.py` 生成：角色限定 driver / oncogene /
tumor_suppressor，按引用数降噪，`#` 注释头记快照日期、来源与过滤口径。

为什么放本地快照：Zenodo 在国内 DNS 被污染，在线下载不可靠；快照让
`--force` 重建完全不联网，同时冻结基准输入、保证可复现。

```bash
# 一次性（VPN 环境）下载 cancermine_collated.tsv 后：
python -m WY.tools.build_pancancer_module --cancermine-file cancermine_collated.tsv --cancermine-update
# 之后（任意网络环境）：
python -m WY.tools.build_pancancer_module --force
```

## Q&A

**Q1：`--force` 报 CancerMine 下载失败 / 超时 / DNS 解析错误，怎么办？**
Zenodo 在国内常被 DNS 污染，不要反复重试。用**能访问 Zenodo 的网络（VPN 或境外主机）**：
1. 打开 Zenodo 上 CancerMine 的记录：<https://doi.org/10.5281/zenodo.1156241>
   （即 CancerMine 论文 Data availability 里给的那个 Zenodo 条目）；
2. 在文件列表里下载 **`cancermine_collated.tsv`**（约 5 MB，CC0）；
3. 回到仓库根目录执行：

   ```bash
   python -m WY.tools.build_pancancer_module --cancermine-file <下载路径>/cancermine_collated.tsv --cancermine-update
   ```

4. 之后 `--force` 只读本地快照，不再访问 Zenodo。

Zenodo 打不开时，[CancerMine 官网](https://cancermine.ai4biomed.org) 的 Download 按钮指向
**同一个 Zenodo 记录**，没有第二官方源；可请能访问的同事代下同一个文件。
**不要**用来路不明的第三方网盘/镜像——无法校验来源，等于污染基准输入。

## 许可

- DailyMed / SPL：公有领域（U.S. Government work）
- Open Targets：CC BY 4.0（署名）
- openFDA：公开域（openFDA terms）
- CancerMine（泛癌集用，见 `WY/tools/build_pancancer_module.py`）：CC0
