# 全局配置文件
import os
from dotenv import load_dotenv

# 自动计算项目根目录（config.py 所在目录，即本地电脑的目录）
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 从 API.env 加载密钥（该文件已被 .gitignore 忽略，不会提交到 git）
load_dotenv(os.path.join(BASE_DIR, "API.env"))

# API密钥
DisGeNET_API_KEY = os.getenv("DisGeNET_API_KEY", "")
BioGRID_API_KEY = os.getenv("BioGRID_API_KEY", "your_api_key")
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
# 疾病靶点路径（DisGeNET 新版 REST API；Header: Authorization=<key>，无 Bearer 前缀）
DisGeNET_URL = "https://api.disgenet.com/api/v1/gda/summary"
# DGIdb v5 GraphQL 端点（药物-基因交互，无需 key）
DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"
# 富集分析路径
gProfiler_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# ChEMBL — 药物靶点主数据源
ChEMBL_API_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
ChEMBL_MAX_ACTIVITIES = 1000
# ======== 备选路径 =========
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
# DGIdb / Open Targets 请求重试次数：这两个 API 单次分析请求量很小（1~2 次），
# 重试过多会在网络不佳时显著拖慢来源更新，2 次足矣（配合退避约多等 3s）
LIGHT_API_MAX_RETRIES = 2
# 富集请求（g:Profiler / Enrichr）重试次数：服务不可达时单次连接可能空等 20s+，
# 重试过多会拖慢整个分析，故单独压低
ENRICHMENT_MAX_RETRIES = 2
# DisGeNET 主动限流：低于官方 TRIAL 配额（每分钟 < 50 次、每天 < 25000 次），
# 发请求前先节流，避免高频请求触发 429 封禁
DISGENET_RATE_LIMIT_PER_MINUTE = 50
DISGENET_RATE_LIMIT_PER_DAY = 25000
# 主动限流等待上限（秒）：每分钟窗口已满时最多等待这么久，超过则拒绝本次请求
DISGENET_RATE_LIMIT_MAX_WAIT = 60
# g:Profiler 熔断标记有效期（秒）：gProfiler 连续失败后写入标记，
# 有效期内直接跳过 g:Profiler 走 Enrichr，避免每次分析都白等连接超时
GPROFILER_CIRCUIT_TTL_SECONDS = 1800
# 缓存 TTL（秒）：参数不变时，PPI / 富集缓存超过该时长视为过期，
# 重新执行同参数分析时重新发请求，以保证拿到最新数据
CACHE_TTL_SECONDS = 6 * 3600
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

# DisGeNET 主动限流状态文件（持久化请求时间戳，跨 Streamlit rerun 生效）
DISGENET_RATE_LIMIT_STATE_FILE = os.path.join(CACHE_DIR, "disgenet_rate_limit.json")

# ======== 日志配置 ========
# 日志目录（WY/Static/log）— 每天一个文件 app_YYYY-MM-DD.log，纯文本 UTF-8
LOG_DIR = os.path.join(BASE_DIR, "Static", "log")
# 日志保留天数：超过该天数的 .log 文件在应用启动时自动删除
LOG_RETENTION_DAYS = 730
# 日志格式：时间 | 级别 | 函数名 | 消息体
# 消息体目标长度 20~120 字（不强制填充/截断，仅作为书写规范）
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

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
# 最小度数过滤器 — 节点度数 ≤ 此值将被过滤，不视为核心靶点，该参数为兜底
MIN_DEGREE_FILTER = 3
# PyVis 物理引擎迭代次数 — 值越大布局越松散，防止节点过度拥挤
PHYSICS_ITERATIONS = 500

# ======== 分子对接配置（AutoDock Vina v1.2.7） ========
# Vina 可执行文件目录（可用环境变量 VINA_DIR 覆盖，便于换机器）
VINA_DIR = os.environ.get("VINA_DIR", r"D:\PythonCode\vina")
VINA_EXE = os.path.join(VINA_DIR, "vina_1.2.7_win.exe")
VINA_SPLIT_EXE = os.path.join(VINA_DIR, "vina_split_1.2.7_win.exe")

# 对接算法参数
EXHAUSTIVENESS = 8          # 搜索深度（越大越充分，耗时越长）
NUM_MODES = 9               # 输出结合模式数量（--num_modes）
DOCKING_CPU = 4             # 并行 CPU 核数（--cpu）
DOCKING_TIMEOUT = 600       # vina subprocess 超时（秒），超时视为失败
DOCKING_SEED = 0            # vina 随机种子（--seed），固定可复现

# 配体 3D 生成参数（方案 B：ChEMBL SMILES + rdkit 生成时使用）
LIGAND_PH = 7.4             # 质子化 pH（rdkit）
LIGAND_3D_SEED = 0          # 构象生成固定种子，保证确定性

# 缓存与归档策略
DOCKING_CACHE_TTL_DAYS = 7  # pdbqt 中间物 / raw 原始文件过期天数，超期自动清理
DOCKING_MAX_VERSIONS = 3    # 同一条件（drug_disease_gene）成果最多保留版本数，旧版入 archive

# ======== 分子对接存储路径 ========
# docking 缓存目录 — 确定性中间物（受体/配体 pdbqt、对接盒），可安全缓存
DOCKING_CACHE_DIR = os.path.join(CACHE_DIR, "docking")
DOCKING_RECEPTOR_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "receptor")
DOCKING_LIGAND_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "ligand")
DOCKING_BOX_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "box")
# raw 目录 — 原始 PDB/SDF，用完即删或随 7 天 TTL 清理
DOCKING_RAW_CACHE_DIR = os.path.join(DOCKING_CACHE_DIR, "raw")

# 对接成果目录 — 每次分析一个独立时间戳子目录，绝对不覆盖
MD_RESULT_DIR = os.path.join(DOWNLOADS_DIR, "MDresult")
# 旧版成果归档目录（超 DOCKING_MAX_VERSIONS 的旧版移入）
ARCHIVE_DIR = os.path.join(MD_RESULT_DIR, "archive")

# ======== 分子对接数据源 ========
# RCSB PDB 下载模板（{pdb_id} 占位符，小写 4 位字母）。
# 说明：受体蛋白 3D 结构只能从 RCSB 下载（pubchempy 只管 PubChem 小分子）；
#       配体 3D SDF 则统一由 pubchempy.get_sdf(cid, record_type="3d") 获取，无需单独 URL。
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
