# 全局配置文件
import logging
import os
from dotenv import load_dotenv

# 平级零依赖叶子:Vina/OBABEL 工具链探测逻辑(本文件导入时调用取值,详见该模块 docstring)。
# 以 WY 包绝对导入(可编辑安装提供可见性,见 selfcheck_architecture)。
from WY import env_discovery

# 项目根目录,即 config.py 所在文件夹
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 从 API.env 读取 API 密钥(该文件不提交到 git)
load_dotenv(os.path.join(BASE_DIR, "API.env"))

# API密钥
DisGeNET_API_KEY = os.getenv("DisGeNET_API_KEY", "")
# ======== 主要路径 =========
# PPI网络路径
# STRING 查询 API 基址不在此写死:按最终生效的 STRING_VERSION 在文件末尾统一派生
# (设置面板可切换版本),保证 PPI 查询与背景互作组始终用同一版本。
# STRING 单次查询基因数上限(约 1000 个)
STRING_QUERY_LIMIT = 1000
# STRING 调用方标识(官方要求,避免被限流)
STRING_CALLER_IDENTITY = "NetworkPharmacologyPlatform"
# 疾病靶点数据源(DisGeNET API),密钥放在 Authorization 请求头中
DisGeNET_URL = "https://api.disgenet.com/api/v1/gda/summary"
# DGIdb 药物-基因交互数据源(无需密钥)
DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"
# ======== API 请求控制(数据源均为国外网站,设置较宽松) ========
API_TIMEOUT = 60
API_MAX_RETRIES = 5
API_RETRY_DELAY = 3
API_REQUEST_DELAY = 0.2
# DGIdb / Open Targets 请求重试次数:每次分析只调 1~2 次,重试 2 次即可。
LIGHT_API_MAX_RETRIES = 2
# ChEMBL 请求超时(秒)与重试次数:靶点详情批量响应大,服务端偶发很慢,放宽以免整个抓取中止
CHEMBL_API_TIMEOUT = 90
CHEMBL_API_MAX_RETRIES = 3
# 富集请求(g:Profiler / Enrichr)尝试次数:连续该次数请求异常(重试机会耗尽)即触发熔断降级,
# 冷却期内直接改用 Enrichr,避免反复苦等超时。
ENRICHMENT_MAX_RETRIES = 3
# DisGeNET 单页请求尝试次数(初次+重试):每次分析仅数页请求,重试 2 次即可。
DISGENET_MAX_RETRIES = 2
# DisGeNET 主动限流:学术(Academic)计划官方 API 配额为 每分钟 50 次 / 每天 25,000 次 /
# 每月 500,000 次。本地限额在官方配额内再留约 10% 余量,避免本地窗口与服务器
# 计数口径差异导致 429。本程序每次分析仅 1 次请求,正常使用远达不到这些上限。
DISGENET_RATE_LIMIT_PER_MINUTE = 45
DISGENET_RATE_LIMIT_PER_DAY = 20000
DISGENET_RATE_LIMIT_PER_MONTH = 450000
# 限流等待上限(秒):超过则放弃本次请求
DISGENET_RATE_LIMIT_MAX_WAIT = 60
# g:Profiler 熔断冷却阶梯(秒):触发熔断从第 1 档开始,之后每次"恢复失败"(冷却到期后的探测失败)
# 推进一档,请求成功(任何有效响应)即清零复位;最后一档为上限。网络故障通常几分钟内恢复,
# 前两档短便于快速恢复,后两档拉长避免频繁撞超时。
GPROFILER_CIRCUIT_COOLDOWN_LADDER = (30, 60, 600, 1800)
# 冷却到期后的探测超时(秒):探测按单次尝试在独立线程限时执行,超时按探测失败继续延长冷却,
# 避免触发探测的分析苦等整个 API_TIMEOUT。
GPROFILER_PROBE_TIMEOUT_SECONDS = 15
# Enrichr 熔断冷却阶梯(秒):兜底源每个请求只尝试 1 次,失败即熔断(重试上移到熔断层,
# 冷却后的探测即重试);首档比 g:Profiler 长,因为每次探测会占用一个模块线程。
ENRICHR_CIRCUIT_COOLDOWN_LADDER = (60, 600, 1800)
# Enrichr 冷却到期后的探测超时(秒),语义与 g:Profiler 相同。
ENRICHR_PROBE_TIMEOUT_SECONDS = 15
# 缓存有效期(秒):PPI/富集/邻近度/RWR 等分析缓存超过该时长即过期,过期只挡复用不删文件。
# 上游源多为月/年级 release,取 7 天减少重复请求;需要立即刷新走管理页"更新"(入口清缓存)。
CACHE_TTL_SECONDS = 7 * 86400
# 数据源 HTTP 响应缓存有效期(小时):四库(ChEMBL/DGIdb/Open Targets/DisGeNET)统一 60 天。
# 60天远短于任一库的发版周期,到期重拉不会错过新数据,又可以把重复请求压到最低;
# 追新走管理页"更新"(入口先清缓存,强制重新联网)。
HTTP_CACHE_TTL_HOURS = 60 * 24
# RCSB 结构搜索的缓存有效期(小时),覆盖全局 HTTP 缓存值:新结构会持续发布,取 4 天折中;
# 结构元数据(GraphQL)条目发布后不可变,沿用全局值即可。
RCSB_SEARCH_CACHE_TTL_HOURS = 96

# ======== 靶点降噪过滤 ========
# 药物/疾病靶点在整合(merge)时做降噪过滤,阈值集中在此配置,过滤逻辑在 core/targets.py。
# 过滤规则或阈值变更后需把 TARGET_FILTER_VERSION 加 1:merged 整合缓存的文件名
# 带该版本号,加 1 后旧缓存自动失效,下次分析按新阈值重建(原始缓存不受影响)。
TARGET_FILTER_VERSION = 3
# ChEMBL 黑名单:非特异性结合/药物代谢蛋白,不计入药物靶点
# (白蛋白 ALB 与主要药物代谢 CYP 酶;CYP46A1 是脑特异胆固醇酶、可能是药效靶点,不列入)
CHEMBL_TARGET_BLACKLIST = (
    "ALB", "CYP1A1", "CYP1A2", "CYP2A6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4",
)
# 设置面板里黑名单的常用候选(白蛋白、主要代谢 CYP、常见结合蛋白);面板也支持手输任意基因符号
CHEMBL_BLACKLIST_OPTIONS = (
    "ALB", "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9",
    "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4", "CYP3A5", "CES1", "GSTP1", "F2", "UGT1A1",
)
# ChEMBL pChEMBL 活性下限:基因的最大 pChEMBL 达到该值才保留;
# 仅有弱定量活性(0 < pChEMBL < 阈值)的基因剔除。
# 完全没有 pChEMBL 定量值的基因默认需被 DGIdb 佐证(通过 DGIDB_MIN_SCORE 过滤)才保留;
# 关闭 CHEMBL_CORROBORATION_REQUIRED 后无定量值的基因无条件保留
# (配合阈值为 0、DGIDB_MIN_SCORE 为 0 即"无降噪")。
# 注意:一刀切按阈值过滤会把靶点砍到个位数(多数活性记录无 pChEMBL 值),勿删 DGIdb 佐证通道
CHEMBL_PCHEMBL_THRESHOLD = 0.5
# 无 pChEMBL 定量值的 ChEMBL 基因是否必须被 DGIdb 佐证才保留(分析参数,可在主页齿轮调整)
CHEMBL_CORROBORATION_REQUIRED = False
# ======== 疾病靶点规模(分析参数,主页齿轮 +- 调节) ========
# Open Targets 整合时取关联分数前 OT_TOP_N 个靶点。原始抓取固定拉满 OT_FETCH_SIZE 条(按分数降序
# 全量落盘),因此调整 Top-N 不需要重新联网;仅当原始缓存行数不足 Top-N 时才自动补抓。
OT_FETCH_SIZE = 1000
OT_TOP_N = 500
# DisGeNET 抓取前 DISGENET_PAGES 页(page_number 分页,学术账号每页固定 100 条,API 硬性规定);
# 整合时按页数×每页条数截断,原始缓存行数不足时自动补抓缺的页。
DISGENET_PAGE_SIZE = 100
DISGENET_PAGES = 5
# DGIdb 交互分数下限,且要求至少 1 条 PMID 文献支撑,两条件同时满足才算有效相互作用
DGIDB_MIN_SCORE = 0.05
# Open Targets 总关联分数下限(保险丝:防少见病的 top200 混入极弱关联靶点);
# 同时只保留蛋白编码(protein_coding)靶点。
OT_MIN_SCORE = 0.1
# ======== 文件存储路径配置 =========

# cache目录 — 储存 API 返回的原始数据
CACHE_DIR = os.path.join(BASE_DIR, "Static", "data", "cache")

# PPI缓存目录 — 储存 STRING API 返回的边表数据
PPI_CACHE_DIR = os.path.join(CACHE_DIR, "ppi")

# 富集缓存目录 — 储存 g:Profiler / Enrichr 返回的富集结果
ENRICHMENT_CACHE_DIR = os.path.join(CACHE_DIR, "enrichment")

# downloads目录 — 储存分析后的成品文件与 zip 压缩包
DOWNLOADS_DIR = os.path.join(BASE_DIR, "Static", "data", "downloads")

# drug目录 — 药物靶点原始缓存（ChEMBL / DGIdb）
DRUG_CACHE_DIR = os.path.join(CACHE_DIR, "drug")

# disease目录 — 疾病靶点原始缓存（Open Targets / DisGeNET）
DISEASE_CACHE_DIR = os.path.join(CACHE_DIR, "disease")

# merged目录 — 去重整合后的药物/疾病靶点缓存
MERGED_CACHE_DIR = os.path.join(CACHE_DIR, "merged")

# merged 整合缓存文件名前缀(带降噪过滤规则版本号;版本升级后旧缓存自动失效重建)
MERGED_DRUG_PREFIX = f"merged_drug_v{TARGET_FILTER_VERSION}_"
MERGED_DISEASE_PREFIX = f"merged_disease_v{TARGET_FILTER_VERSION}_"

# DisGeNET 限流状态文件(保存请求时间戳,页面刷新后仍生效)
DISGENET_RATE_LIMIT_STATE_FILE = os.path.join(CACHE_DIR, "disgenet_rate_limit.json")

# 数据源 HTTP 响应缓存(DGIdb/Open Targets/DisGeNET)的 sqlite 落盘位置
HTTP_CACHE_FILE = os.path.join(CACHE_DIR, "http_cache.sqlite")

# CLI 结果目录:取 API.env 里的 CLI_RESULTS_DIR(配一次永久生效,与 VINA_DIR 同类的本机路径,不入库);
# 未配置时为空串,cli.py 兜底落到 ./results。供 cli.py 的 --out 三级解析作第二优先级。
CLI_RESULTS_DIR = os.getenv("CLI_RESULTS_DIR", "").strip()

# ======== 日志配置 ========
# 日志目录:每天一个日志文件 app_YYYY-MM-DD.log
LOG_DIR = os.path.join(BASE_DIR, "Static", "log")
# 日志保留天数(天):超期日志在应用启动时自动删除
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
LOG_RETENTION_DAYS = 365
# 日志格式:时间 | 级别 | 函数名 | 消息
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# ======== 物种映射表 =========
# 目前收录人、小鼠、大鼠,需要其他物种可自行添加
SPECIES_MAP = {
    "Human": {
        "ncbi_taxon": "9606",       # 给 STRING / BioGRID 用
        "gprofiler": "hsapiens",    # 给 g:Profiler 用
    },
    "Mouse": {
        "ncbi_taxon": "10090",
        "gprofiler": "mmusculus",
    },
    "Rat": {
        "ncbi_taxon": "10116",
        "gprofiler": "rnorvegicus",
    }
}
# ======== 分析参数 =========
# PPI 网络置信度阈值(STRING):0.4 为中等置信度(多数文献取值),0.7 为高置信度
PPI_CONFIDENCE = 0.4
# 富集分析显著性阈值
ENRICHMENT_PVALUE = 0.055
# 富集分析的条目类型:GO 生物过程 / 分子功能 / 细胞组分,以及 KEGG 通路
ENRICHMENT_SOURCES = ("GO:BP", "GO:MF", "GO:CC", "KEGG")
# 每类富集展示条数上限
ENRICHMENT_TOP_PER_SOURCE = 10
# 通路展示数量上限
PATHWAY_DISPLAY_NUM = 30
# CytoHubba 式核心靶点:按 MCC(Maximal Clique Centrality)排序取前 N 个画核心靶点子网络图
MCC_TOP_N = 10

# ======== Hub 基因筛选(排名制 top-N + 小网络保护) ========
# 排名依据:"degree"(对齐文献按度数取前 N 的惯例)或 "mcc"(cytoHubba 默认算法)
HUB_METHOD = "degree"
# 排名制取前 N 个作为 Hub(文献惯例为 cytoHubba top 10)
HUB_TOP_N = 10
# 小网络保护系数:节点数 < HUB_SMALL_NET_FACTOR*HUB_TOP_N 时 top-N 覆盖面过大,
# 回退为阈值法(度 > 全部节点中位数)
HUB_SMALL_NET_FACTOR = 2
# 网络图布局迭代次数:值越大节点分布越散,不易挤在一起
PHYSICS_ITERATIONS = 500

# ======== 蛋白模块检测（Louvain + 模块富集） ========
# Louvain 分辨率与随机种子(种子固定保证结果可复现)
LOUVAIN_RESOLUTION = 1.0
LOUVAIN_SEED = 0
# 参与富集的模块规模下限与数量上限:基因太少的模块不做富集
LOUVAIN_MIN_MODULE_SIZE = 5
LOUVAIN_MAX_ENRICH_MODULES = 10
# 每个模块页面展示的富集条目上限(完整结果见 CSV 下载)
LOUVAIN_TOP_TERMS_PER_MODULE = 5
# 若最大模块占网络比例超过该值,页面会提示调整分辨率
LOUVAIN_MAX_MODULE_RATIO = 0.6

# ======== 排序富集（GSEA prerank，仅智人） ========
# Enrichr 基因集库:首次使用时自动下载保存到本地,之后不再重复下载
# (删除 gmt 缓存目录即可重新下载)。默认只取 Hallmark + KEGG,需要时可自行添加。
GSEA_ENRICHR_LIBS = ("MSigDB_Hallmark_2020", "KEGG_2021_Human")
# GMT 快照目录
GSEA_GMT_DIR = os.path.join(CACHE_DIR, "gmt")
# 基因集规模范围:基因数在此范围外的通路不参与分析
GSEA_SET_MIN_SIZE = 15
GSEA_SET_MAX_SIZE = 500
# 置换次数:结果已缓存,想要更稳的 p 值可调高后重跑
GSEA_PERMUTATION_NUM = 1000
# 随机种子(保证结果可复现)
GSEA_SEED = 0
# 排序基因数下限,低于该值不计算(结果无统计意义)
GSEA_MIN_RANKED_GENES = 30
# 页面展示条数(按 |NES| 从大到小)
GSEA_TOP_DISPLAY = 20
# 页面显著性筛选阈值(按惯例 FDR<0.25;仅用于筛选展示,完整结果见 CSV)
GSEA_FDR_SCREEN = 0.25

# ======== 背景互作组(STRING 全网络本地缓存) ========
# 查询 API 基址与下载路径都由本版本派生,保证 PPI 查询与背景网络计算用同一版本数据;
# 可在设置面板切换(写 ui_settings.json),新版本发布后无需改代码。
STRING_VERSION = "v12.0"
# 预设版本下拉:官方下载站在售的版本(2026-09 实测:v11.0/v11.5/v12.0/v12.5,v10.5 及更早已下线);
# 更新的版本直接手输即可(经下载站探活后接受)
STRING_KNOWN_VERSIONS = ("v12.0", "v12.5")


def normalized_string_version(value) -> str | None:
    """校验并归一 STRING 版本号:合法返回 'v<数字>[.<数字>][字母后缀]' 小写形式(如 v12.0、
    v12.0a、v10;漏 v 前缀自动补),非法返回 None。设置面板保存与 ui_settings.json 加载共用,
    保证写盘与读回的格式一致。允许字母后缀:官方历史上有补丁版命名(如 v10.0a)。"""
    import re

    s = str(value).strip().lower()
    if s and not s.startswith("v"):
        s = f"v{s}"
    return s if re.fullmatch(r"v\d+(\.\d+)?[a-z]?", s) else None


# 官方下载站
STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download"
# 官方静态备用站(与主站同数据,支持 Range 分段)
STRING_DOWNLOAD_STATIC_BASE = "https://stringdb-static.org/download"
# 国内镜像:2026-10 起下载直链已全面停用(仅保留网页门户,页面文件链接也指向官方站),
# 留在链尾作哨兵——官方若恢复镜像直链,把它在 STRING_DOWNLOAD_SOURCES 里挪回首 位即可。
STRING_DOWNLOAD_MIRROR_BASE = "https://cn.string-db.org/download"
# 下载源尝试顺序:官方站优先(国内无 STRING 数据镜像,清华 TUNA 等仅镜像 Bioconductor
# 的 R 包,不含数据文件);官方站单连接实测仅 ~10KB/s,多段并行可近似线性提速。
INTERACTOME_DOWNLOAD_SOURCES = (
    STRING_DOWNLOAD_BASE, STRING_DOWNLOAD_STATIC_BASE, STRING_DOWNLOAD_MIRROR_BASE)
# 全网络大文件下载并发分段数(仅对支持 HTTP Range 的源生效);不支持 Range 或分段
# 失败自动回退单流下载并提醒用户(慢,约 10 分钟以上)。
INTERACTOME_DOWNLOAD_CONCURRENCY = 6
# 背景互作组缓存目录
INTERACTOME_CACHE_DIR = os.path.join(CACHE_DIR, "interactome")
# 置信度换算:PPI_CONFIDENCE 是 0~1 小数,全网络打分是 0~1000 整数。
# 派生值统一在本文件末尾(UI 设置覆盖之后)计算,保证与最终生效的 PPI_CONFIDENCE 一致。
# 背景网络目前只构建人类;小鼠/大鼠需要时再放开(数据量约为 3 倍)
INTERACTOME_ENABLED_SPECIES = ("Human",)
# 背景网络"陈旧提醒"阈值(天):超过该天数只在页面提示检查新版本,不会自动重新下载
# (同版本文件内容不变;真正的更新是等 STRING 出新版本后改 STRING_VERSION 并重建)。
INTERACTOME_STALE_DAYS = 365

# ======== 网络邻近度检验(依赖背景互作组) ========
# 随机置换次数:结果已缓存,想要更稳的结果可调高
PROXIMITY_N_PERM = 200
# 随机种子（结果可复现）
PROXIMITY_SEED = 0
# 度匹配分箱宽度:连接数相近的节点归入同一箱、随机替换(按论文 Algorithm 1 实现)。
PROXIMITY_DEGREE_BIN_SIZE = 100
# 药物/疾病靶点映射到背景网络的比例低于该值时拒绝计算
PROXIMITY_MIN_COVERAGE = 0.3
# 邻近度结果缓存目录(只在计算成功时保存)
PROXIMITY_CACHE_DIR = os.path.join(CACHE_DIR, "proximity")

# ======== RWR 靶点优先级排序(依赖背景互作组) ========
# 重启概率 α:越大游走越远、排序越偏全局,越小越贴近已知靶点;0.7 为常用取值。
RWR_ALPHA = 0.7
# 收敛阈值与最大迭代次数
RWR_TOL = 1e-8
RWR_MAX_ITER = 200
# 页面展示的非种子基因数量(已知靶点单独全量展示,不占此数)
RWR_TOP_N = 50
# RWR 扩展靶点数:分析时取 RWR 得分最高的前 N 个非种子基因与交集靶点取并集
# 作为查询集(0 = 关闭扩展)。默认 30,与文献交叉靶点规模对齐。
RWR_EXPAND_TOP_N = 30
# 缓存中保留的非种子排名行数(已知靶点全量另存,不占此名额)。
RWR_CACHE_TOP_N = 500
# 度匹配校正检验的置换次数与种子:检验 Top 排名是否只是高连接度带来的假象。
RWR_NULL_PERMS = 200
RWR_NULL_SEED = 0
# RWR 结果缓存目录(只在计算成功时保存)
RWR_CACHE_DIR = os.path.join(CACHE_DIR, "rwr")

# ======== 分子对接配置（AutoDock Vina v1.2.7） ========
# Vina 程序目录自动探测(逻辑在 env_discovery.py,导入时求值):优先环境变量 VINA_DIR
# (可写入 API.env),其次 PATH;换电脑只需设置其中之一。都找不到时 VINA_DIR 为空串,
# 运行对接时会报错提示设置 VINA_DIR。
VINA_DIR = env_discovery.find_vina_dir()
VINA_EXE = env_discovery.find_exe(VINA_DIR, "vina", exclude="vina_split")
VINA_SPLIT_EXE = env_discovery.find_exe(VINA_DIR, "vina_split")  # 缺失时 docking 自动改用内置拆分
# 受体转换以子进程调用 obabel 命令行(不依赖其 Python 绑定);都找不到时 OBABEL_EXE 为空串,运行转换时报错提示。
OBABEL_DIR = env_discovery.find_obabel_dir()
OBABEL_EXE = env_discovery.find_exe(OBABEL_DIR, "obabel")

# 对接算法参数
EXHAUSTIVENESS = 8          # 搜索深度（越大越充分，耗时越长）
NUM_MODES = 9               # 输出结合模式数量（--num_modes）
DOCKING_CPU = 4             # 并行 CPU 核数（--cpu）,0 表示自动用全部核心
DOCKING_TIMEOUT = 600       # 单次对接超时(秒),超时视为失败
OBABEL_TIMEOUT = 120        # 单次 obabel 转换超时(秒),超时视为失败
DOCKING_SEED = 0            # vina 随机种子（--seed）
PUBCHEM_BUSY_RETRIES = 2    # PubChem 503 繁忙自动重试次数(瞬时抖动常见,总尝试=次数+1)
PUBCHEM_BUSY_BACKOFF = 2.0  # 重试退避基数(秒),逐次加倍(2s、4s)

# 缓存与归档策略
DOCKING_CACHE_TTL_DAYS = 7  # 加工后的中间文件(受体/配体 pdbqt 等)保留天数,超期自动清理
DOCKING_RAW_TTL_DAYS = 1    # 原始/过程文件(SDF/PDB 等)保留天数,过期清理后可重新生成
DOCKING_MAX_VERSIONS = 3    # 同一药物-疾病组合的成果最多保留版本数,旧版直接删除

# ======== 分子对接存储路径 ========
# 对接缓存目录:受体/配体 pdbqt、对接盒等中间文件
DOCKING_CACHE_DIR = os.path.join(CACHE_DIR, "docking")
DOCKING_RECEPTOR_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "receptor")
DOCKING_LIGAND_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "ligand")
# 原始文件目录:原始 PDB/SDF,到期自动清理
DOCKING_RAW_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "raw")

# 对接成果目录:每次分析单独建时间戳文件夹,不覆盖旧结果
MD_RESULT_DIR = os.path.join(DOWNLOADS_DIR, "MDresult")
# 收藏目录:点击 ☆ 收藏的成果复制到此,不会被自动清理
MD_STAR_DIR = os.path.join(DOWNLOADS_DIR, "star")
# GEO+ML 筛选成果目录(页面):与 MDresult 同级的 Web 端产物区;
# CLI 默认仍落 results/mlscreen(方案 §6),页面历史列表同时扫两处
GEO_RESULT_DIR = os.path.join(DOWNLOADS_DIR, "GEOresult")

# ======== 分子对接数据源 ========
# 受体蛋白 PDB 结构下载地址({pdb_id} 为 4 位 PDB 编号);配体 SDF 由 PubChem 获取。
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# ======== 分子对接:AlphaFold 预测结构兜底 ========
# 总开关:还需在页面上另行开启(页面开关默认关闭)才使用 AlphaFold 兜底
ALPHAFOLD_ENABLED = True
# AlphaFold 模型版本:官方升级时只需改这里,旧缓存(文件名含版本)自然失效。
ALPHAFOLD_MODEL_VERSION = "v6"
# 结构文件下载地址({accession} 为 UniProt 主条目号,{version} 为上面的模型版本)
ALPHAFOLD_FILE_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_{version}.pdb"
# 基因名→UniProt 号的查询渠道,按顺序依次尝试,某个不可用自动换下一个。
UNIPROT_LOOKUP_SOURCES = ("ebi", "uniprot")
EBI_PROTEINS_API_URL = "https://www.ebi.ac.uk/proteins/api/proteins"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
# 基因名→UniProt 号查询结果的缓存文件(长期有效)
UNIPROT_LOOKUP_CACHE = os.path.join(DOCKING_CACHE_DIR, "uniprot_gene_map.json")
# AlphaFold 结构没有共晶配体,对接盒初始取整个蛋白的外包围盒;任一边超过该值(埃)
# 时自动向中心收拢。仅为可用起点,页面会提示人工确认口袋并微调。
DOCKING_BOX_MAX_SIZE_AF = 40.0
# RCSB 检索失败重试次数;RCSB 元数据下载并发数
RCSB_SEARCH_MAX_RETRIES = 3
RCSB_DATA_MAX_WORKERS = 4
# 对接口袋残基判定半径(埃):配体周围该距离内的受体残基视为结合位点残基,文献常用 4 Å。
# 属于分析参数,可在主页齿轮(Analysis settings)调整,持久化到 ui_settings.json。
POCKET_CUTOFF = 4.0

# ======== 各数据源最小请求间隔(秒) ========
# 只对需要连续快速请求的站点做限速,未列出的数据源请求量小,不做额外限制。
HOST_MIN_REQUEST_INTERVAL = {
    "chembl": 0.2,        # ChEMBL 批量分页:保持原 API_REQUEST_DELAY 的礼貌间隔
    "rcsb_search": 0.3,   # RCSB 批量检索:约 3 次/秒,配合 4 路并发仍属礼貌范围
    "rcsb_data": 0.25,    # RCSB 元数据下载
    "ebi_proteins": 1.1,  # EBI 官方建议 ≤1 次/秒（UniProt 基因查询）
    "alphafold": 1.1,     # AlphaFold 文件下载（EBI 同域）
    "uniprot": 1.1,       # UniProt REST 兜底
}

# ======== 页面设置覆盖(主页面板齿轮保存的运行时设置) ========
# 齿轮面板修改的分析参数会写入 Static/data/cache/ui_settings.json,
# 启动时在此套用,保证重启后设置仍生效。只允许覆盖白名单内的键。
UI_SETTINGS_FILE = os.path.join(CACHE_DIR, "ui_settings.json")
# 搜索历史持久化文件:各页面搜索框的历史按搜索框 key 分组存放,重启仍生效
SEARCH_HISTORY_FILE = os.path.join(CACHE_DIR, "search_history.json")
# 每个搜索框最多保留的历史条数(最新的在前)
SEARCH_HISTORY_MAX = 6
UI_SETTINGS_KEYS = (
    "CHEMBL_PCHEMBL_THRESHOLD", "DGIDB_MIN_SCORE", "OT_MIN_SCORE",
    "ENRICHMENT_PVALUE", "LOUVAIN_RESOLUTION", "POCKET_CUTOFF",
    "CHEMBL_CORROBORATION_REQUIRED", "OT_TOP_N", "DISGENET_PAGES",
    "PPI_CONFIDENCE",  # STRING 置信度阈值:PPI 边筛选与背景互作组共用
    "LOG_LEVEL",  # Logs 页的 DEBUG 开关:持久化后重启按此级别初始化日志
    "STRING_VERSION",  # STRING 数据版本:查询 API/下载路径/各缓存键随版本派生或隔离
    "UI_FONT_SCALE",  # 前端界面等比缩放(1.0=100%):纯显示偏好,不影响任何分析计算
)
# UI_FONT_SCALE 允许范围:加载设置文件时校验,坏值按内置默认运行
_UI_FONT_SCALE_RANGE = (0.8, 1.5)
# 界面等比缩放默认 100%;实际生效值由 ui_settings.json 覆盖(下方 apply_ui_settings_overrides)
UI_FONT_SCALE = 1.0
# LOG_LEVEL 允许的取值:加载设置文件时校验,坏值按内置默认运行
_LOG_LEVEL_NAMES = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")


def apply_ui_settings_overrides() -> None:
    """读取 ui_settings.json,把白名单内的运行时设置覆盖到本模块常量。

    设置文件损坏或键值非法时按本文件内置默认值运行,不阻断启动;
    selfcheck 可改写 UI_SETTINGS_FILE 后重调本函数验证加载逻辑。
    """
    if not os.path.exists(UI_SETTINGS_FILE):
        return
    try:
        import json
        with open(UI_SETTINGS_FILE, "r", encoding="utf-8") as _f:
            _ui_overrides = json.load(_f)
        for _k, _v in _ui_overrides.items():
            if _k == "CHEMBL_TARGET_BLACKLIST":
                # 黑名单是字符串列表,单独校验后套用
                if isinstance(_v, list) and all(isinstance(_i, str) for _i in _v):
                    globals()[_k] = tuple(_v)
            elif _k == "STRING_VERSION":
                # 版本号经统一校验归一,坏值按内置默认运行;必须先于通用数值/布尔分支,
                # 防止手改设置文件写入数字/布尔被原样赋值后,末尾派生 URL 时崩溃
                _sv = normalized_string_version(_v)
                if _sv:
                    globals()[_k] = _sv
            elif _k == "UI_FONT_SCALE":
                # 缩放比例限范围,越界/坏值按内置默认 1.0 运行;先于通用数值分支
                if isinstance(_v, (int, float)) and not isinstance(_v, bool) \
                        and _UI_FONT_SCALE_RANGE[0] <= _v <= _UI_FONT_SCALE_RANGE[1]:
                    globals()[_k] = float(_v)
            elif _k in UI_SETTINGS_KEYS and isinstance(_v, bool):
                globals()[_k] = _v
            elif _k in UI_SETTINGS_KEYS and isinstance(_v, (int, float)) \
                    and not isinstance(_v, bool):
                globals()[_k] = _v
            elif _k == "LOG_LEVEL" and isinstance(_v, str) \
                    and _v.strip().upper() in _LOG_LEVEL_NAMES:
                globals()[_k] = _v.strip().upper()
    except Exception:
        # 设置文件损坏时不阻断启动,按本文件内置值运行
        logging.getLogger(__name__).warning(
            "ui_settings.json 读取失败,按内置默认设置运行", exc_info=True)


# import 时套用一次:对使用方而言行为与原先"import 即生效"完全一致
apply_ui_settings_overrides()

# PPI 置信度可能被上面的设置文件覆盖,而背景互作组的置信度阈值在上方按内置值派生:
# 覆盖生效后必须重算,保证两处阈值始终一致(运行时齿轮面板修改亦需同样同步)
INTERACTOME_CONFIDENCE = int(PPI_CONFIDENCE * 1000)

# STRING 版本可能被上面的设置文件覆盖:查询 API 基址按最终生效版本派生
# (STRING 版本化子域名把 v12.0 写作 version-12-0),与下载路径/背景互作组保持同版本;
# 运行时齿轮面板切换版本亦需同步重算此值
STRING_BASE_URL = f"https://version-{STRING_VERSION[1:].replace('.', '-')}.string-db.org/api/json"


# ======== GEO / ML 筛选(GEO与ML共识筛选模块,实施方案 §5) ========
# E-utilities 官方限流:无 key 3 次/秒,注册免费 NCBI_API_KEY 后 10 次/秒;
# key 放 API.env(可选,缺失自动按无 key 模式限速,功能不受影响)
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_SEARCH_CACHE_TTL_HOURS = 24      # 检索结果短期缓存(每天都有新数据集入库)
GEO_SERIES_CACHE_SUBDIR = "geo_series"  # series matrix 永久缓存子目录(GEO 归档不可变)
GEO_DOWNLOAD_TIMEOUT = 600           # 锁接管兜底窗口基准(超时×重试+60s);下载已改分段续传,无整体超时
GEO_DOWNLOAD_RETRIES = 3             # 下载失败重试次数(间隔见下),gzip 损坏视为失败重试
GEO_DOWNLOAD_RETRY_DELAY = 5         # 重试间隔(秒)
GEO_DOWNLOAD_CONCURRENCY = 4         # 大文件分段并发连接数:跨境对 NCBI 单连接常被压到个位数 KB/s
                                     # (实测 8.5MB 爬 16 分钟),Range 分段多连接聚合提速;E-utilities 的
                                     # 3 req/s 限流只约束检索接口不约束文件服务器,4 连接在容忍范围内
GEO_DOWNLOAD_STALL_TIMEOUT = 90      # 单连接断流判定(秒):连续无任何字节到达即断开该段;重试按段续传
ML_SCREENING_VERSION = 1             # 算法/权重变更时 +1,旧结果自动失效(仿 TARGET_FILTER_VERSION)
ML_MIN_SAMPLES_PER_GROUP = 3         # 每组样本低于此数直接报错,不降级硬跑
ML_SEED = 0                          # 全模块固定随机种子,可复现
ML_BOOTSTRAP_ROUNDS = 30             # bootstrap 稳定性重采样轮次(2026-10 实测定档:500基因×30轮 ML 段约3.5分钟,页面/CLI 统一)
ML_CV_FOLDS = 5                      # 交叉验证折数(M2)
ML_SHORTLIST_TOP_K = 20              # shortlist 截断(M2)
ML_DEG_PREFILTER_TOP_N = 500         # 进 ML 前按 FDR 截断的候选基因数(性能前提;500 与 30 轮为实测耗时定档,页面/CLI 统一)
# NCBI 限流注册:按是否配 key 取对应间隔(方案:不新建限流器,并入现有节流表)
HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] = 0.10 if NCBI_API_KEY else 0.34

ML_CONSENSUS_RANK_W = 0.5             # 共识打分权重:平均秩分(方案:权重写死 config,杜绝调出来)
ML_CONSENSUS_FREQ_W = 0.5             # 共识打分权重:bootstrap 选中频率
ML_RESULTS_KEEP_PER_GSE = 6           # 结果保留语义(§10 增补):同 tag 只留最近 6 份,落盘后自动清
                                      # 理更旧的;跨数据集总量永不自动清(实验记录),治理走结果管理页


# ======== API.env 可视化设置(齿轮面板:数据源板块的密钥 + 分子对接板块的本机工具路径) ========
# 键名与 API.env 文件内一致;落盘由 core.management.save_api_env 负责(set_key 保留注释与其他行),
# 落盘后经 apply_api_env 回填运行时,免重启生效。密钥类键不回显当前值,输入新值即覆盖。
API_ENV_KEYS = ("DisGeNET_API_KEY", "NCBI_API_KEY", "VINA_DIR", "OBABEL_DIR", "CLI_RESULTS_DIR")


def apply_api_env(changes: dict) -> list[str]:
    """把用户提交的 API.env 参数回填到运行时 config 常量,返回实际生效的键名列表。

    只接受 API_ENV_KEYS 白名单内的键;派生值随源同步重算:NCBI_API_KEY 变更联动
    E-UTILS 节流间隔(与本文件前文注册逻辑保持同一变换),VINA/OBABEL 目录变更时
    重探对应可执行文件(目录为空时 find_exe 安全返回空串;VINA 目录变更时
    VINA_SPLIT_EXE 一并重探)。本函数只改本进程状态,不做任何落盘。"""
    global DisGeNET_API_KEY, NCBI_API_KEY, VINA_DIR, VINA_EXE, VINA_SPLIT_EXE
    global OBABEL_DIR, OBABEL_EXE, CLI_RESULTS_DIR
    applied: list[str] = []
    for _k, _v in changes.items():
        if _k not in API_ENV_KEYS:
            continue
        _v = str(_v).strip()
        if _k == "DisGeNET_API_KEY":
            DisGeNET_API_KEY = _v
        elif _k == "NCBI_API_KEY":
            NCBI_API_KEY = _v
            HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] = 0.10 if _v else 0.34
        elif _k == "VINA_DIR":
            VINA_DIR = _v
            VINA_EXE = env_discovery.find_exe(_v, "vina", exclude="vina_split")
            VINA_SPLIT_EXE = env_discovery.find_exe(_v, "vina_split")
        elif _k == "OBABEL_DIR":
            OBABEL_DIR = _v
            OBABEL_EXE = env_discovery.find_exe(_v, "obabel")
        elif _k == "CLI_RESULTS_DIR":
            CLI_RESULTS_DIR = _v
        else:
            continue
        applied.append(_k)
    return applied
