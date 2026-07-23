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
# ======== API请求控制(鉴于都是国外网站所以设置的条件更宽松) =========
API_TIMEOUT = 60
API_MAX_RETRIES = 5
API_RETRY_DELAY = 3
API_REQUEST_DELAY = 1
# ======== 文件存储路径配置（先不写，讨论之后再写） =========
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
#生信的我不懂，后面找你了解
