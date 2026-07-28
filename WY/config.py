# 全局配置文件
# API密钥
DisGeNET_API_KEY = "your_api_key"
BioGRID_API_KEY = "your_api_key"
# ======== 主要路径 =========
# 药物靶点路径
STITCH_URL = "https://stitch.embl.de/api/json/getInteractions"
# PPI网络路径
STRING_URL = "https://string-db.org/api/json/network"
# 疾病靶点路径
DisGeNET_URL = "https://www.disgenet.org/api/gda/disease"
# 富集分析路径
gProfiler_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# ======== 备选路径 =========
# 药物靶点替补路径
PubChem_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
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
API_REQUEST_DELAY = 1
# ======== 文件存储路径配置 =========
import os
# 自动计算项目根目录（config.py 所在目录本地电脑的目录）
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# cache目录 — 储存 API 返回的原始数据
CACHE_DIR = os.path.join(BASE_DIR, "Static", "data", "cache")

# downloads目录 — 储存分析后的成品文件与 zip 压缩包
DOWNLOADS_DIR = os.path.join(BASE_DIR, "Static", "data", "downloads")

# ======== 物种映射表 =========
#先写了人类，小鼠，大鼠，后面想加自己加
SPECIES_MAP = {
    "Human": {
        "ncbi_taxon": "9606",       # 给 STRING / STITCH / BioGRID 用
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
# 药物-靶点连接阈值（STITCH，范围 0~1），≥0.7 为可信
DRUG_TARGET_SCORE = 0.7
# 通路展示数量上限
PATHWAY_DISPLAY_NUM = 30
# 最小度数过滤器 — 节点度数 ≤ 此值将被过滤，不视为核心靶点，该参数为兜底，正常应该使用中位数法
MIN_DEGREE_FILTER = 3
# PyVis 物理引擎迭代次数 — 值越大布局越松散，防止节点过度拥挤
PHYSICS_ITERATIONS = 500
