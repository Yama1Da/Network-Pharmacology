# 全局配置文件
import os
from dotenv import load_dotenv

# 项目根目录,即 config.py 所在文件夹
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 从 API.env 读取 API 密钥(该文件不提交到 git)
load_dotenv(os.path.join(BASE_DIR, "API.env"))

# API密钥
DisGeNET_API_KEY = os.getenv("DisGeNET_API_KEY", "")
# ======== 主要路径 =========
# PPI网络路径
# 固定版本的 STRING 地址,保证分析结果可复现
STRING_BASE_URL = "https://version-12-0.string-db.org/api/json"
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
# 富集请求(g:Profiler / Enrichr)重试次数:服务不可达时等待时间长,故压低。
ENRICHMENT_MAX_RETRIES = 2
# DisGeNET 主动限流:学术(Academic)计划官方 API 配额为 每分钟 50 次 / 每天 25,000 次 /
# 每月 500,000 次。本地限额在官方配额内再留约 10% 余量,避免本地窗口与服务器
# 计数口径差异导致 429。本程序每次分析仅 1 次请求,正常使用远达不到这些上限。
DISGENET_RATE_LIMIT_PER_MINUTE = 45
DISGENET_RATE_LIMIT_PER_DAY = 20000
DISGENET_RATE_LIMIT_PER_MONTH = 450000
# 限流等待上限(秒):超过则放弃本次请求
DISGENET_RATE_LIMIT_MAX_WAIT = 60
# g:Profiler 连续失败后的熔断时长(秒):期间直接改用 Enrichr,避免反复等超时。
GPROFILER_CIRCUIT_TTL_SECONDS = 1800
# 缓存有效期(秒):PPI/富集缓存超过该时长即过期,下次分析重新请求。
CACHE_TTL_SECONDS = 6 * 3600
# 数据源 HTTP 响应缓存有效期(小时):DGIdb/Open Targets/DisGeNET 共用,与 ChEMBL 库默认一致
HTTP_CACHE_TTL_HOURS = 24

# ======== 靶点降噪过滤 ========
# 药物/疾病靶点在整合(merge)时做降噪过滤,阈值集中在此配置,过滤逻辑在 core/targets.py。
# 过滤规则或阈值变更后需把 TARGET_FILTER_VERSION 加 1:merged 整合缓存的文件名
# 带该版本号,加 1 后旧缓存自动失效,下次分析按新阈值重建(原始缓存不受影响)。
TARGET_FILTER_VERSION = 1
# ChEMBL 黑名单:非特异性结合/药物代谢蛋白,不计入药物靶点
# (白蛋白 ALB 与主要药物代谢 CYP 酶;CYP46A1 是脑特异胆固醇酶、可能是药效靶点,不列入)
CHEMBL_TARGET_BLACKLIST = (
    "ALB", "CYP1A2", "CYP2A6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4",
)
# ChEMBL pChEMBL 活性下限:基因的最大 pChEMBL 达到该值才保留;
# 完全没有 pChEMBL 定量值的基因需被 DGIdb 佐证(通过 DGIDB_MIN_SCORE 过滤)才保留;
# 仅有弱定量活性(0 < pChEMBL < 阈值)的基因剔除。
# 注意:一刀切按阈值过滤会把靶点砍到个位数(多数活性记录无 pChEMBL 值),勿删 DGIdb 佐证通道
CHEMBL_PCHEMBL_THRESHOLD = 4.8
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
# PPI 网络置信度阈值(STRING):0.7 为高置信度,0.4 为中置信度
PPI_CONFIDENCE = 0.7
# 富集分析显著性阈值
ENRICHMENT_PVALUE = 0.05
# 富集分析的条目类型:GO 生物过程 / 分子功能 / 细胞组分,以及 KEGG 通路
ENRICHMENT_SOURCES = ("GO:BP", "GO:MF", "GO:CC", "KEGG")
# 每类富集展示条数上限
ENRICHMENT_TOP_PER_SOURCE = 10
# 通路展示数量上限
PATHWAY_DISPLAY_NUM = 30
# 最小度数过滤:连接数不超过该值的节点不作为核心靶点(兜底过滤)
MIN_DEGREE_FILTER = 3
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

# ======== 背景互作组(STRING v12.0 全网络本地缓存) ========
# 版本需与上面的 STRING API 一致,保证 PPI 查询与背景网络计算用同一版本数据。
STRING_VERSION = "v12.0"
# 官方下载站(人类全网络文件约 83MB)
STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download"
# 国内镜像:优先使用,不可用时自动改用官方站。
STRING_DOWNLOAD_MIRROR_BASE = "https://cn.string-db.org/download"
# 背景互作组缓存目录
INTERACTOME_CACHE_DIR = os.path.join(CACHE_DIR, "interactome")
# 置信度换算:上面 PPI_CONFIDENCE 是 0~1 小数,全网络打分是 0~1000 整数,由这里换算。
INTERACTOME_CONFIDENCE = int(PPI_CONFIDENCE * 1000)  # 0.7 -> 700
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
# 缓存中保留的非种子排名行数(已知靶点全量另存,不占此名额)。
RWR_CACHE_TOP_N = 500
# 度匹配校正检验的置换次数与种子:检验 Top 排名是否只是高连接度带来的假象。
RWR_NULL_PERMS = 200
RWR_NULL_SEED = 0
# RWR 结果缓存目录(只在计算成功时保存)
RWR_CACHE_DIR = os.path.join(CACHE_DIR, "rwr")

# ======== 分子对接配置（AutoDock Vina v1.2.7） ========
# Vina 程序目录自动探测:优先环境变量 VINA_DIR(可写入 API.env),其次 PATH;换电脑只需设置其中之一。
# 都找不到时 VINA_DIR 为空串,运行对接时会报错提示设置 VINA_DIR。
def _find_exe(vina_dir: str, prefix: str, exclude: str = "") -> str:
    """目录里按文件名前缀找可执行文件(不限定平台与版本文件名)。

    多个匹配时取字典序最后一个;exclude 用于主程序排除 vina_split;找不到返回空串。
    """
    if not vina_dir or not os.path.isdir(vina_dir):
        return ""
    names = sorted(
        f for f in os.listdir(vina_dir)
        if f.lower().startswith(prefix)
        and not (exclude and f.lower().startswith(exclude))
        and os.path.isfile(os.path.join(vina_dir, f))
    )
    return os.path.join(vina_dir, names[-1]) if names else ""


def _find_vina_dir() -> str:
    """探测 Vina 程序所在目录(跨平台):VINA_DIR 优先,其次 PATH。

    VINA_DIR 或 PATH 里填目录即可(也可填可执行文件全路径)。
    """
    env_dir = os.environ.get("VINA_DIR")
    if env_dir:
        if os.path.isfile(env_dir):                                   # 填的是可执行文件全路径
            return os.path.dirname(env_dir)
        if _find_exe(env_dir, "vina", exclude="vina_split"):          # 填的是目录
            return env_dir

    # vina 所在目录已加入 PATH(按 vina* 前缀匹配,兼容各平台发布文件名)
    for path_dir in os.environ.get("PATH", "").split(os.pathsep):
        if _find_exe(path_dir.strip().strip('"'), "vina", exclude="vina_split"):
            return path_dir.strip().strip('"')

    return ""  # 找不到时,运行对接会报错提示设置 VINA_DIR


VINA_DIR = _find_vina_dir()
VINA_EXE = _find_exe(VINA_DIR, "vina", exclude="vina_split")
VINA_SPLIT_EXE = _find_exe(VINA_DIR, "vina_split")  # 缺失时 docking 自动改用内置拆分

# 对接算法参数
EXHAUSTIVENESS = 8          # 搜索深度（越大越充分，耗时越长）
NUM_MODES = 9               # 输出结合模式数量（--num_modes）
DOCKING_CPU = 4             # 并行 CPU 核数（--cpu）
DOCKING_TIMEOUT = 600       # 单次对接超时(秒),超时视为失败
DOCKING_SEED = 0            # vina 随机种子（--seed），固定可复现

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
