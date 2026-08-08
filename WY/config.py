# 全局配置文件
# API密钥
DisGeNET_API_KEY = "your_api_key"
BioGRID_API_KEY = "your_api_key"
# ======== 主要路径 =========
# PPI网络路径
# 版本化基础地址：裸 string-db.org 会随 STRING 版本更新而变化，版本化地址保证分析结果可复现
STRING_BASE_URL = "https://version-12-0.string-db.org/api/json"
# 非版本化地址（可手动切换，一般不推荐用于正式分析）
STRING_URL = "https://string-db.org/api/json/network"
# STRING 单次查询基因数量上限（API 限制约 1000 个）
STRING_QUERY_LIMIT = 1000
# STRING 调用方标识（STRING 官方要求标注应用身份，避免被限流）
STRING_CALLER_IDENTITY = "NetworkPharmacologyPlatform"
# 疾病靶点路径
DisGeNET_URL = "https://www.disgenet.org/api/gda/disease"
# 富集分析路径
gProfiler_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# ChEMBL — 药物靶点主数据源
ChEMBL_API_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
ChEMBL_MAX_ACTIVITIES = 1000
# ======== 备选路径 =========
# PPI网络替补路径
BioGRID_URL = "https://webservice.thebiogrid.org/interactions"
# 疾病靶点替补路径
OpenTargetsPlatform_URL = "https://api.platform.opentargets.org/api/v4/graphql"
# 富集分析替补路径
ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_QUERY_URL = "https://maayanlab.cloud/Enrichr/enrich"
# ======== 其他路径 =========
# RxNorm — 药物名称转 RxCUI，使用方式：requests.get(RxNorm_URL, params=RxNorm_PARAMS(drug_name))
RxNorm_URL = "https://rxnav.nlm.nih.gov/REST/rxcui.json"
def RxNorm_PARAMS(drug_name: str) -> dict:
    return {"idtype": "name", "name": drug_name}
# ======== API请求控制(鉴于都是国外网站所以设置的条件更宽松) =========
API_TIMEOUT = 60
API_MAX_RETRIES = 5
API_RETRY_DELAY = 3
API_REQUEST_DELAY = 0.2
# 富集请求（g:Profiler / Enrichr）重试次数：服务不可达时单次连接可能空等 20s+，
# 重试过多会拖慢整个分析，故单独压低
ENRICHMENT_MAX_RETRIES = 2
# g:Profiler 熔断标记有效期（秒）：gProfiler 连续失败后写入标记，
# 有效期内直接跳过 g:Profiler 走 Enrichr，避免每次分析都白等连接超时
GPROFILER_CIRCUIT_TTL_SECONDS = 1800
# 缓存 TTL（秒）：参数不变时，PPI / 富集缓存超过该时长视为过期，
# 重新执行同参数分析时重新发请求，以保证拿到最新数据
CACHE_TTL_SECONDS = 6 * 3600
# ======== 文件存储路径配置 =========
import os
# 自动计算项目根目录（config.py 所在目录本地电脑的目录）
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# cache目录 — 储存 API 返回的原始数据
CACHE_DIR = os.path.join(BASE_DIR, "Static", "data", "cache")

# PPI缓存目录 — 储存 STRING API 返回的边表数据
PPI_CACHE_DIR = os.path.join(CACHE_DIR, "ppi")

# 富集缓存目录 — 储存 g:Profiler / Enrichr 返回的富集结果
ENRICHMENT_CACHE_DIR = os.path.join(CACHE_DIR, "enrichment")

# downloads目录 — 储存分析后的成品文件与 zip 压缩包
DOWNLOADS_DIR = os.path.join(BASE_DIR, "Static", "data", "downloads")

# ======== 物种映射表 =========
#先写了人类，小鼠，大鼠，后面想加自己加
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
# PPI 网络置信度阈值（STRING），≥0.7 为高置信度，≥0.4为中置信度
PPI_CONFIDENCE = 0.7
# 富集分析显著性阈值
ENRICHMENT_PVALUE = 0.05
# 富集分析的 sources（g:Profiler 命名）：
# GO:BP 生物过程 / GO:MF 分子功能 / GO:CC 细胞组分 / KEGG 通路
ENRICHMENT_SOURCES = ("GO:BP", "GO:MF", "GO:CC", "KEGG")
# 每个 source 的富集展示条数上限
ENRICHMENT_TOP_PER_SOURCE = 10
# 药物-靶点连接阈值（Open Targets），≥0.7 为可信
DRUG_TARGET_SCORE = 0.7
# 通路展示数量上限
PATHWAY_DISPLAY_NUM = 30
# 最小度数过滤器 — 节点度数 ≤ 此值将被过滤，不视为核心靶点，该参数为兜底，正常应该使用中位数法
MIN_DEGREE_FILTER = 3
# PyVis 物理引擎迭代次数 — 值越大布局越松散，防止节点过度拥挤
PHYSICS_ITERATIONS = 500
