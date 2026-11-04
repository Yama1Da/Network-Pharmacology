"""PPI 网络分析:从 STRING 获取蛋白互作边表,构建网络并计算拓扑指标,绘图并保存结果。

绘图渲染器拆在 ppi_plots(核心靶点图、pyvis 交互图与各风格静态图)。
"""

import hashlib
import logging
import os

import pandas as pd
import requests

from .caching import (clean_expired_cache_files, is_cache_fresh, open_csv_cache,
                      read_cache_key, save_csv_cache, save_json_file)
from .http_client import request_with_retry
from .text_utils import genes_from_dataframe, safe_name
# 图渲染器已拆分到 ppi_plots(纯渲染函数,图数据由调用方传入);
# 再导入以维持 WY.core.ppi 既有对外接口(analysis/tests 从本模块导入),build_main_figures 亦经此调用。
from .ppi_plots import (_ball_sprite, _edge_evidence_map, render_core_target_plot,
                        render_cytoscape_style_plot, render_pyvis_html,
                        render_static_plot, render_string_style_plot)


def _string_post(endpoint: str, params: dict) -> requests.Response | None:
    """向 STRING 数据库发送请求,失败时返回 None。
    提交的基因在 STRING 中一个都识别不出来(404),不算错误,按空结果处理,不重试。"""
    from WY import config

    url = f"{config.STRING_BASE_URL}/{endpoint}"
    resp = request_with_retry(
        url, data=params, max_retries=config.API_MAX_RETRIES,
        timeout=config.API_TIMEOUT, label=f"STRING {endpoint}",
        no_retry_status=(404,),
    )
    if resp is not None and resp.status_code == 404:
        logging.warning("STRING %s: 404 (no identifiers resolved)", endpoint)
        return None
    return resp


def fetch_string_ppi(
    genes: list[str],
    species: str = "9606",
    score: float = 0.7,
    network_type: str = "full",
    label: str | None = None,
) -> pd.DataFrame:
    """从 STRING 数据库获取基因列表的蛋白互作边表:先把基因符号映射为 STRING 蛋白 ID(不认识的基因被剔除),再按置信度阈值 score(0~1)筛选互作边。
    network_type 默认 "full" 为功能关联,"physical" 仅保留有物理结合证据的互作;species 为物种编号,人类是 "9606";基因列表自动去重排序,最多取 1000 个。
    相同参数的结果存为本地缓存,有效期内且基因列表指纹与文件头记录一致时直接复用,否则重新请求覆盖;请求失败、基因全部不被识别或没有达到阈值的边时,返回空表。"""
    from WY import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        return pd.DataFrame()
    if len(genes) > config.STRING_QUERY_LIMIT:
        genes = genes[:config.STRING_QUERY_LIMIT]

    # 基因列表短指纹写入缓存文件头,命中时校验:缓存文件名只有参数不含基因列表,
    # 上游变化(设置调整后交集/RWR 扩展更新等)会在 TTL 内让同名缓存命中旧查询集的边表,
    # 校验指纹不符立即视为未命中重算覆盖(旧格式无指纹行的文件按不匹配处理)。
    gene_md5 = hashlib.md5("\r".join(genes).encode("utf-8")).hexdigest()
    # 缓存键含 STRING 版本:切版本后同名缓存按键不符处理,不会在 TTL 内命中旧版本边表
    # (旧格式无版本前缀的文件同样按未命中,重算覆盖一次)
    cache_key = f"{config.STRING_VERSION}|{gene_md5[:12]}"

    # 缓存文件名由来源标识、网络类型、置信度阈值拼成;network_type 必须参与命名,
    # 否则 full/physical 切换后 TTL 内会命中另一参数组合的旧边表。
    # 相同参数组合固定同名,重新生成时直接覆盖旧文件。
    if label:
        label_part = safe_name(label)
        name = f"StringPPI_{label_part}_{network_type}_{score:g}.csv"
        if label_part == "manual":
            # 手动基因模式的基因列表无法从名字区分,附加基因列表短指纹作缓存标识,
            # 避免不同基因列表混用同一份缓存
            name = f"StringPPI_manual_{network_type}_{score:g}_{gene_md5[:8]}.csv"
    else:
        name = f"StringPPI_{network_type}_{score:g}.csv"
    cache_file = os.path.join(config.PPI_CACHE_DIR, name)

    # TTL 过期文件在命中判定里只会被跳过、不会被删除,这里顺手清理避免缓存目录只增不减
    clean_expired_cache_files(config.PPI_CACHE_DIR)

    # 命中需同时满足 TTL 未过期与文件头指纹一致(与富集缓存同一判定顺序):
    # 先验键再读文件,指纹不符或旧格式无键行按未命中处理,避免白读文件与误导性命中日志
    cached = None
    if is_cache_fresh(cache_file) and read_cache_key(cache_file) == cache_key:
        cached = open_csv_cache(cache_file, fresh=True, label=f"STRING PPI for {len(genes)} genes")
    if cached is not None:
        return cached

    id_resp = _string_post("get_string_ids", {
        "identifiers": "\r".join(genes),
        "species": species,
        "caller_identity": config.STRING_CALLER_IDENTITY,
    })
    if id_resp is None:
        return pd.DataFrame()
    try:
        id_rows = id_resp.json()
    except Exception:
        logging.exception("STRING get_string_ids response parse failed")
        return pd.DataFrame()

    string_ids = [r.get("stringId") for r in id_rows if r.get("stringId")]
    if not string_ids:
        logging.warning("STRING: none of %d genes were recognized", len(genes))
        return pd.DataFrame()

    net_resp = _string_post("network", {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "required_score": int(score * 1000),
        "network_type": network_type,
        "caller_identity": config.STRING_CALLER_IDENTITY,
    })
    if net_resp is None:
        return pd.DataFrame()
    try:
        rows = net_resp.json()
    except Exception:
        logging.exception("STRING network response parse failed")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    if df.empty:
        logging.info("STRING: no edges above score %s for %d genes", score, len(genes))
        return df

    gene_set = {g.upper() for g in genes}
    df = df[
        df["preferredName_A"].str.upper().isin(gene_set)
        & df["preferredName_B"].str.upper().isin(gene_set)
    ].rename(columns={"preferredName_A": "protein_A", "preferredName_B": "protein_B"})

    keep = ["protein_A", "protein_B", "score", "nscore", "fscore", "pscore",
            "ascore", "escore", "dscore", "tscore"]
    df = df[[c for c in keep if c in df.columns]]

    save_csv_cache(df, cache_file, header_lines=[cache_key])
    logging.info("STRING PPI for %d genes: %d edges cached", len(genes), len(df))
    return df


def build_ppi_graph(edges_df: pd.DataFrame, all_genes: list[str]) -> "nx.Graph":
    """用 STRING 边表构建蛋白互作网络图。
    所有输入基因都作为节点(没有互作关系的基因成为孤立节点);同一对基因的重复边只保留置信度最高的一条。"""
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(all_genes)
    if edges_df.empty:
        return G

    edges = edges_df[edges_df["protein_A"] != edges_df["protein_B"]].copy()
    edges["_key"] = edges.apply(
        lambda r: tuple(sorted((str(r["protein_A"]), str(r["protein_B"])))), axis=1
    )
    edges = edges.sort_values("score", ascending=False).drop_duplicates(subset="_key")
    for _, r in edges.iterrows():
        s = float(r["score"])
        G.add_edge(r["protein_A"], r["protein_B"], score=s, weight=s / 1000.0)
    return G


def compute_ppi_metrics(graph: "nx.Graph") -> pd.DataFrame:
    """计算每个节点的拓扑指标:度(degree)、介数中心性、接近中心性;孤立节点各项为 0,结果按度数从高到低排列。"""
    import networkx as nx

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph)
    closeness = nx.closeness_centrality(graph)

    df = pd.DataFrame({
        "gene": list(graph.nodes()),
        "degree": [degree[n] for n in graph.nodes()],
        "betweenness": [betweenness.get(n, 0.0) for n in graph.nodes()],
        "closeness": [closeness.get(n, 0.0) for n in graph.nodes()],
    })
    return df.sort_values("degree", ascending=False).reset_index(drop=True)


def mark_hub_genes(metrics_df: pd.DataFrame,
                   top_n: int = 10, method: str = "degree",
                   mcc_map: dict[str, float] | None = None,
                   small_net_factor: int = 2) -> pd.DataFrame:
    """标记 Hub 基因(排名制 top-N + 小网络保护),在原表基础上新增三列:is_hub(True 为 Hub 基因)、hub_rank(按排名依据的位次,1 起;并列按度数、基因名稳定排序)、hub_rule(本次采用的规则,排名制为 "top{N}",小网络回退为 "threshold")。
    规则:节点数 >= small_net_factor*top_n 时,按 method 指定的拓扑量(degree 或 mcc)取前 top_n 个——与文献按 cytoHubba top-10 筛选核心靶点的惯例一致;网络过小时 top-N 覆盖面过大(超过 1/small_net_factor),回退为阈值法:连接数高于全部节点中位数(建网时孤立节点已被隐藏,中位数 >= 1,故等价于要求至少两条互作)。
    method="mcc" 但 mcc_map 缺失/为空时自动回退 degree 排名(MCC 计算失败或图过度稀疏的兜底)。返回表按排名依据降序重排。"""
    df = metrics_df.copy()
    n = len(df)
    if n == 0:
        df["is_hub"] = []
        df["hub_rank"] = []
        df["hub_rule"] = []
        return df

    eff_method = method if method in ("degree", "mcc") else "degree"
    use_rank = n >= max(2, int(small_net_factor)) * int(top_n)

    if use_rank and eff_method == "mcc":
        if not mcc_map or not (df["gene"].map(mcc_map).fillna(0.0) > 0).any():
            eff_method = "degree"  # MCC 缺失/全零(稀疏图)时退回度排名

    if use_rank:
        if eff_method == "mcc":
            df["mcc"] = df["gene"].map(mcc_map).fillna(0.0)
            df = df.sort_values(["mcc", "degree", "gene"],
                                ascending=[False, False, True], kind="mergesort")
        else:
            df = df.sort_values(["degree", "gene"],
                                ascending=[False, True], kind="mergesort")
        df["hub_rank"] = range(1, n + 1)
        df["is_hub"] = df["hub_rank"] <= int(top_n)
        df["hub_rule"] = f"top{int(top_n)}"
    else:
        threshold = df["degree"].median()
        df["is_hub"] = df["degree"] > threshold
        df = df.sort_values(["degree", "gene"],
                            ascending=[False, True], kind="mergesort")
        df["hub_rank"] = range(1, n + 1)
        df["hub_rule"] = "threshold"
    return df.reset_index(drop=True)


# MCC 的最大团枚举数量上限:超过按失败处理(由调用方跳过核心靶点图),
# 防止病态稠密网络把 maximal clique 枚举拖成无限等待
_MCC_CLIQUE_LIMIT = 500_000


def compute_mcc_scores(graph: "nx.Graph") -> dict[str, float]:
    """计算每个节点的最大团中心性 MCC(Maximal Clique Centrality,CytoHubba 默认算法):MCC(v) = Σ_{u∈N(v)} (C(v,u) - 1)!,
    其中 C(v,u) 为同时包含节点 v 与邻居 u 的最大团(maximal clique)的大小;没有三角形的节点退化为 MCC = 度数,孤立节点为 0。
    实现上对全图做一次最大团枚举(nx.find_cliques),对每对相邻节点记录同时包含它们的所有最大团里的最大规模,再按公式聚合;图的最大团数量超过 _MCC_CLIQUE_LIMIT 时抛 RuntimeError。"""
    import networkx as nx
    from math import factorial

    # pair_clique[v][u] = 目前见过的同时包含 v、u 的最大团规模
    pair_clique: dict[str, dict[str, int]] = {n: {} for n in graph.nodes()}
    n_cliques = 0
    for clique in nx.find_cliques(graph):
        n_cliques += 1
        if n_cliques > _MCC_CLIQUE_LIMIT:
            raise RuntimeError(
                f"maximal cliques exceed {_MCC_CLIQUE_LIMIT}, MCC computation aborted")
        k = len(clique)
        if k < 2:
            continue  # 孤立点构成的长度 1 团没有 (v,u) 节点对
        for v in clique:
            cand = pair_clique[v]
            for u in clique:
                if u != v and cand.get(u, 0) < k:
                    cand[u] = k
    return {
        v: float(sum(factorial(c - 1) for c in nbrs.values()))
        for v, nbrs in pair_clique.items()
    }


def pick_top_mcc_genes(metrics_df: pd.DataFrame, mcc_map: dict[str, float],
                       top_n: int = 10) -> list[str]:
    """按 MCC 从高到低取前 top_n 个基因,作为核心靶点;MCC 并列时按度数、再按基因名稳定排序。
    MCC 为 0(无任何互作)的节点不参与入选——全无互作时返回空列表,调用方据此跳过核心靶点图。"""
    df = metrics_df[["gene", "degree"]].copy()
    df["mcc"] = df["gene"].map(mcc_map).fillna(0.0)
    df = df[df["mcc"] > 0]
    df = df.sort_values(["mcc", "degree", "gene"],
                        ascending=[False, False, True], kind="mergesort")
    return df.head(top_n)["gene"].tolist()


def build_main_figures(graph, metrics_df, edges_df, hub_method: str, mcc_map: dict):
    """渲染 pyvis 交互图、静态图、STRING/Cytoscape 风格图与核心靶点图
    (pyplot 全局状态非线程安全,调用方须在单线程内完成本函数)。

    MCC 得分缺失时在此补算,并入 metrics_df(随节点 CSV 落盘);核心靶点图与两张风格图
    失败只跳过自身。返回 (html_content, static_fig, mcc_map, hub_df, core_fig,
    core_genes, string_fig, cyto_fig)。
    """
    from WY import config

    html_content = render_pyvis_html(graph, metrics_df)
    static_fig = render_static_plot(graph, metrics_df)
    string_fig = cyto_fig = None
    try:
        string_fig = render_string_style_plot(graph, metrics_df, edges_df)
    except Exception:
        logging.exception("STRING-style figure failed — skipped")
    try:
        cyto_fig = render_cytoscape_style_plot(graph, metrics_df)
    except Exception:
        logging.exception("Cytoscape-style figure failed — skipped")
    if not mcc_map:
        try:
            mcc_map = compute_mcc_scores(graph)
        except Exception:
            logging.exception("MCC computation failed — core target figure skipped")
    hub_df = metrics_df[metrics_df["is_hub"]]  # 默认值:MCC 失败时 hub 表不带 mcc 列
    if mcc_map:
        metrics_df["mcc"] = metrics_df["gene"].map(mcc_map)
        hub_df = metrics_df[metrics_df["is_hub"]]  # 重取让 Hub 表带上 mcc 列
    core_fig = None
    core_genes: list[str] = []
    if mcc_map:
        core_genes = pick_top_mcc_genes(metrics_df, mcc_map, config.MCC_TOP_N)
        if core_genes:
            try:
                core_fig = render_core_target_plot(graph, core_genes, mcc_map)
            except Exception:
                logging.exception("Core target figure failed — skipped")
    return (html_content, static_fig, mcc_map, hub_df, core_fig, core_genes,
            string_fig, cyto_fig)


def base_disease_of(name: str) -> str:
    """剥掉历史实验变体后缀得到干净疾病名(rwrNN,不区分大小写)。

    旧分析曾把 RWR 参数变体标签拼进疾病名输入框以隔离结果目录
    (如 "gastric cancer rwr72"),下拉选疾病/检索 GEO 时应使用剥掉后缀的
    干净名,变体差异只体现在结果目录与其中的排序 CSV 上。
    """
    import re
    return re.sub(r"\s*rwr\s*\d+\s*$", "", str(name), flags=re.IGNORECASE).strip()


def list_analyzed_pairs() -> list[dict]:
    """列出本机已跑过的 PPI 分析(即 App 页面分析过的药病对),目录名倒序。

    数据源:downloads 下各 <prefix>_data/analysis_meta.json 的 drug/disease 键
    (2026-10 起的分析必带,精确);更早的旧结果无此二键,退回从目录名前缀反解
    —— prefix 由 f"{drug}_{disease}" 经 re.sub 把非单词字符压成 _ 而来,反解为
    best-effort:首个 _ 前作药名、其余作病名、_ 还原为空格。"manual" 目录是
    手动基因输入模式,无药病对,跳过。

    base_disease 为剥掉 rwrNN 变体后缀的干净疾病名(ML 页疾病下拉的选项来源);
    disease 保留原样(含变体后缀,供区分实验批次)。string_version 供排序 CSV
    发现处展示。供 ML 筛选页构建疾病下拉与排序文件上下文匹配(方案 §M4)。
    返回 [{drug, disease, base_disease, string_version, analyzed_at,
           source: "meta"|"prefix", directory}]。
    """
    import json

    from WY import config
    out: list[dict] = []
    dl = config.DOWNLOADS_DIR
    if not os.path.isdir(dl):
        return out
    for d in sorted(os.listdir(dl), reverse=True):
        dpath = os.path.join(dl, d)
        if not os.path.isdir(dpath) or not d.endswith("_data"):
            continue
        prefix = d[:-len("_data")]
        if prefix.startswith("manual"):
            # 手动基因桶:目录名是基因列表哈希,没有真实药病语境,不进 ML 页的
            # 疾病发现与排序候选(曾因表单残留值冒充 metformin/T2DM 被发现过)
            continue
        meta = {}
        meta_path = os.path.join(dpath, "analysis_meta.json")
        if os.path.isfile(meta_path):
            try:
                with open(meta_path, encoding="utf-8") as f:
                    meta = json.load(f)
            except Exception:
                meta = {}
        analyzed_at = str(meta.get("analyzed_at", ""))
        if meta.get("disease"):
            disease = str(meta["disease"])
            out.append({"drug": str(meta.get("drug", "")),
                        "disease": disease,
                        "base_disease": base_disease_of(disease),
                        "string_version": str(meta.get("string_version", "")),
                        "analyzed_at": analyzed_at, "source": "meta",
                        "directory": dpath})
        else:
            head, sep, tail = prefix.partition("_")
            if not sep or not tail:
                continue
            disease = tail.replace("_", " ")
            out.append({"drug": head, "disease": disease,
                        "base_disease": base_disease_of(disease),
                        "string_version": str(meta.get("string_version", "")),
                        "analyzed_at": analyzed_at, "source": "prefix",
                        "directory": dpath})
    return out


def list_analyzed_diseases() -> list[str]:
    """去重后的已分析疾病清单(剥 rwrNN 后缀的干净名,最近分析在前)。

    供 ML 筛选页的疾病下拉直接使用(方案 §M4);大小写不敏感去重,
    保留首次出现的写法。"""
    seen: set[str] = set()
    out: list[str] = []
    for rec in list_analyzed_pairs():
        d = rec["base_disease"].strip()
        if d and d.lower() not in seen:
            seen.add(d.lower())
            out.append(d)
    return out


def list_analyzed_drugs() -> list[str]:
    """去重后的已分析药物清单(最近分析在前)。

    供 ML 筛选页排序文件的药物过滤下拉使用;大小写不敏感去重,
    保留首次出现的写法。"""
    seen: set[str] = set()
    out: list[str] = []
    for rec in list_analyzed_pairs():
        d = rec["drug"].strip()
        if d and d.lower() not in seen:
            seen.add(d.lower())
            out.append(d)
    return out


def list_ranked_csv_candidates() -> list[dict]:
    """列出全部已分析药病对下的排序 CSV(*_query_genes.csv,gene/source/rwr_score
    排序表),供 ML 筛选页勾选进入 shortlist 交集(方案 §M4:按元数据匹配,
    不用文件名猜)。

    每项 {label, path, drug, base_disease};label = 目录前缀 · 分析时间 ·
    STRING 版本。nodes/edges 是网络结构表而非排序,不在此列。"""
    out: list[dict] = []
    for rec in list_analyzed_pairs():
        try:
            dir_files = os.listdir(rec["directory"])
        except OSError:
            continue
        for f in sorted(dir_files):
            if not f.endswith("_query_genes.csv"):
                continue
            label = " · ".join(x for x in (
                os.path.basename(rec["directory"])[:-len("_data")],
                rec["analyzed_at"][:19],
                f"STRING {rec['string_version']}" if rec["string_version"] else "",
            ) if x)
            out.append({"label": label,
                        "path": os.path.join(rec["directory"], f),
                        "drug": rec["drug"],
                        "base_disease": rec["base_disease"]})
    return out


def save_ppi_results(
    graph: "nx.Graph",
    static_fig: "matplotlib.figure.Figure",
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    html_content: str,
    prefix: str,
    enrichment_df: pd.DataFrame | None = None,
    enrichment_bar_fig: "matplotlib.figure.Figure | None" = None,
    enrichment_dot_fig: "matplotlib.figure.Figure | None" = None,
    core_fig: "matplotlib.figure.Figure | None" = None,
    target_pathway_fig: "matplotlib.figure.Figure | None" = None,
    string_fig: "matplotlib.figure.Figure | None" = None,
    cyto_fig: "matplotlib.figure.Figure | None" = None,
    drug: str | None = None,
    disease: str | None = None,
    genes: list | None = None,
) -> dict:
    """把一次 PPI 分析的全部结果保存到结果目录下的 <prefix>_data 文件夹(prefix 为文件名前缀,如 aspirin_cancer),返回各文件的保存路径(含目录路径)。
    保存前会删除同 prefix 的旧文件夹(含历史结果),同一药物+疾病组合只保留最新一次分析,不会越积越多;节点表含 is_hub 列(算出 MCC 时含 mcc 列),其余为边表、交互式网络、GraphML、静态图(PNG 300 dpi + SVG 矢量)、STRING/Cytoscape 风格网络图及各富集/核心靶点/通路网络图表。
    传入富集结果表/图时才保存对应的富集文件,不传则跳过。
    传入 drug/disease 时写入 analysis_meta.json,供下游按元数据匹配(方案 §M4,不用文件名猜)。"""
    import glob
    import shutil
    import networkx as nx
    from WY import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    for old in glob.glob(os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data*")):
        if os.path.isdir(old):
            shutil.rmtree(old)
    subdir = os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data")
    os.makedirs(subdir, exist_ok=True)

    # 结果自描述:记录分析时的数据版本与阈值,任何一份旧结果都能说明自己用的
    # STRING 版本(切版本后对比/复现的依据);结果目录每次分析整体重建,本文件随之刷新
    from datetime import datetime

    meta = {
        "string_version": config.STRING_VERSION,
        "interactome_confidence": int(config.INTERACTOME_CONFIDENCE),
        "analyzed_at": datetime.now().isoformat(timespec="seconds"),
    }
    if drug and str(drug).strip():
        meta["drug"] = str(drug).strip()
    if disease and str(disease).strip():
        meta["disease"] = str(disease).strip()
    # 手动基因模式的目录名后缀是基因列表哈希,人读不出内容——把分析基因清单
    # 原样写进 meta,作为桶的唯一人类可读入口(names 模式已有 drug/disease 可识别,不重复记录)
    if genes and str(prefix).startswith("manual"):
        meta["genes"] = sorted({str(g).upper() for g in genes})
    save_json_file(os.path.join(subdir, "analysis_meta.json"), meta,
                   indent=2, sort_keys=True)

    paths = {
        "nodes": os.path.join(subdir, f"{prefix}_nodes.csv"),
        "edges": os.path.join(subdir, f"{prefix}_edges.csv"),
        "html": os.path.join(subdir, f"{prefix}_network.html"),
        "graphml": os.path.join(subdir, f"{prefix}_graph.graphml"),
    }
    nodes_df.to_csv(paths["nodes"], index=False, encoding="utf-8-sig")
    edges_df.to_csv(paths["edges"], index=False, encoding="utf-8-sig")
    with open(paths["html"], "w", encoding="utf-8") as f:
        f.write(html_content)
    nx.write_graphml(graph, paths["graphml"])
    paths["png"], paths["svg"] = _save_fig_pair(static_fig, subdir, prefix, "network")

    if enrichment_df is not None and not enrichment_df.empty:
        paths["enrichment_csv"] = os.path.join(subdir, f"{prefix}_enrichment.csv")
        enrichment_df.to_csv(paths["enrichment_csv"], index=False, encoding="utf-8-sig")
    if enrichment_bar_fig is not None:
        paths["enrichment_bar_png"], paths["enrichment_bar_svg"] = _save_fig_pair(
            enrichment_bar_fig, subdir, prefix, "enrichment_bar")
    if enrichment_dot_fig is not None:
        paths["enrichment_dot_png"], paths["enrichment_dot_svg"] = _save_fig_pair(
            enrichment_dot_fig, subdir, prefix, "enrichment_dot")
    if core_fig is not None:
        paths["core_png"], paths["core_svg"] = _save_fig_pair(
            core_fig, subdir, prefix, "core_targets")
    if target_pathway_fig is not None:
        paths["dtp_png"], paths["dtp_svg"] = _save_fig_pair(
            target_pathway_fig, subdir, prefix, "enrichment_network")
    if string_fig is not None:
        paths["string_png"], paths["string_svg"] = _save_fig_pair(
            string_fig, subdir, prefix, "network_string")
    if cyto_fig is not None:
        paths["cyto_png"], paths["cyto_svg"] = _save_fig_pair(
            cyto_fig, subdir, prefix, "network_cytoscape")

    paths["directory"] = subdir
    return paths


def _save_fig_pair(fig: "matplotlib.figure.Figure", subdir: str,
                   prefix: str, stem: str) -> tuple[str, str]:
    """把图保存为 <prefix>_<stem>.png + .svg 一对,返回 (png 路径, svg 路径)。
    PNG 固定 300 dpi,满足期刊对位图的分辨率要求;SVG 为矢量,缩放无损。"""
    png = os.path.join(subdir, f"{prefix}_{stem}.png")
    svg = os.path.join(subdir, f"{prefix}_{stem}.svg")
    fig.savefig(png, format="png", dpi=300, bbox_inches="tight")
    fig.savefig(svg, format="svg", bbox_inches="tight")
    return png, svg


def _save_optional_csvs(directory: str, label: str,
                        specs: tuple[tuple[str, str, pd.DataFrame | None], ...]) -> dict:
    """把若干结果表追加保存到既有分析目录(<prefix>_data),specs 每项为 (结果键, 文件名, 表格);
    约定:空表不落盘,避免产生无意义的空文件。返回实际写出的文件路径,调用方会并入结果字典,打包下载时自动包含。"""
    os.makedirs(directory, exist_ok=True)
    paths: dict = {}
    for key, name, df in specs:
        if df is not None and not df.empty:
            path = os.path.join(directory, name)
            df.to_csv(path, index=False, encoding="utf-8-sig")
            paths[key] = path
    if paths:
        logging.info("%s saved to %s (%d files)", label, directory, len(paths))
    return paths


def save_module_results(
    directory: str,
    prefix: str,
    summary_df: pd.DataFrame,
    enrichment_df: pd.DataFrame | None = None,
) -> dict:
    """把蛋白模块检测结果追加保存到既有分析目录(<prefix>_data):包括模块汇总表({prefix}_modules.csv)和
    模块富集长表({prefix}_module_enrichment.csv)。约定与 _save_optional_csvs 相同:空表不落盘。"""
    return _save_optional_csvs(directory, "Module results", (
        ("modules_csv", f"{prefix}_modules.csv", summary_df),
        ("module_enrichment_csv", f"{prefix}_module_enrichment.csv", enrichment_df),
    ))


def save_module_computation(directory: str, prefix: str, computed: dict) -> dict:
    """模块计算结果落盘并附上保存路径(自动流程与交互区手动重跑共用);
    空结果不保存。保存到分析目录后,打包下载（ZIP）会自动包含模块文件。"""
    mod_saved = save_module_results(directory, prefix,
                                    computed["summary_df"], computed["long_df"])
    return {**computed, "saved_paths": mod_saved}


def read_ppi_dataset(dir_name: str) -> tuple[set[str], set[str]]:
    """读取指定的 PPI 结果目录,返回 HUB 基因集合与富集验证基因集合;
    对应文件缺失时返回空集合。历史版本产出的 nodes 表可能缺 is_hub/gene 列
    (downloads 目录跨版本留存),缺列或读不了时按空集合处理,不让调用方崩溃。"""
    from WY import config

    base = os.path.join(config.DOWNLOADS_DIR, dir_name)
    prefix = dir_name[:-5] if dir_name.endswith("_data") else dir_name
    hub: set[str] = set()
    nodes_path = os.path.join(base, f"{prefix}_nodes.csv")
    if os.path.exists(nodes_path):
        try:
            nodes = pd.read_csv(nodes_path, encoding="utf-8-sig", comment="#")
            if {"is_hub", "gene"}.issubset(nodes.columns):
                hub = set(nodes.loc[nodes["is_hub"].fillna(False), "gene"]
                          .dropna().astype(str).str.upper())
            else:
                logging.warning("nodes file lacks is_hub/gene columns, hub list "
                                "skipped: %s", nodes_path)
        except Exception:
            logging.exception("Failed to read nodes file: %s", nodes_path)
    enr: set[str] = set()
    enr_path = os.path.join(base, f"{prefix}_enrichment.csv")
    if os.path.exists(enr_path):
        try:
            enr = genes_from_dataframe(pd.read_csv(enr_path, encoding="utf-8-sig", comment="#"))
        except Exception:
            logging.exception("Failed to read enrichment file: %s", enr_path)
    return hub, enr


def save_gsea_results(
    directory: str,
    prefix: str,
    drug_df: pd.DataFrame | None = None,
    disease_df: pd.DataFrame | None = None,
) -> dict:
    """把 GSEA prerank 结果追加保存到既有分析目录(<prefix>_data):包括药物靶点排秩结果({prefix}_gsea_drug.csv)
    和疾病靶点排秩结果({prefix}_gsea_disease.csv),以及各侧的 NES 条形图({prefix}_gsea_{side}.png/.svg)。
    条形图只画与页面展示一致的 FDR 筛后 top-N,该侧筛后为空则不画也不落盘——无图即无文件,
    ZIP 与下载按钮自动缺省;CSV 始终是完整计算结果。约定同前:空表不落盘。"""
    from WY import config
    from .enrichment import render_gsea_plot

    paths = _save_optional_csvs(directory, "GSEA results", (
        ("gsea_drug_csv", f"{prefix}_gsea_drug.csv", drug_df),
        ("gsea_disease_csv", f"{prefix}_gsea_disease.csv", disease_df),
    ))
    for side, df in (("drug", drug_df), ("disease", disease_df)):
        if df is None or df.empty or "fdr_q_value" not in df.columns:
            continue
        screened = df[df["fdr_q_value"] < config.GSEA_FDR_SCREEN]
        fig = render_gsea_plot(screened) if not screened.empty else None
        if fig is None:
            continue
        png, svg = _save_fig_pair(fig, directory, prefix, f"gsea_{side}")
        paths[f"gsea_{side}_png"] = png
        paths[f"gsea_{side}_svg"] = svg
    return paths


def save_rwr_results(
    directory: str,
    prefix: str,
    disease_df: pd.DataFrame | None = None,
    drug_df: pd.DataFrame | None = None,
) -> dict:
    """把 RWR 靶点优先级排序结果追加保存到既有分析目录(<prefix>_data):包括以疾病靶点为种子的排序表
    ({prefix}_rwr_disease_seeds.csv)和以药物靶点为种子的排序表({prefix}_rwr_drug_seeds.csv)。约定同前:空表不落盘。"""
    return _save_optional_csvs(directory, "RWR results", (
        ("rwr_disease_csv", f"{prefix}_rwr_disease_seeds.csv", disease_df),
        ("rwr_drug_csv", f"{prefix}_rwr_drug_seeds.csv", drug_df),
    ))


def save_query_genes(directory: str, prefix: str,
                     expansion_df: pd.DataFrame | None) -> dict:
    """保存分析查询集清单({prefix}_query_genes.csv),gene/source/rwr_score 三列。RWR 扩展启用时记录每个基因来自
    种子(seed)还是扩展(expanded)及对应得分,保证富集查询集可审计、可复现;未启用扩展(expansion_df 为 None)不落盘。"""
    if expansion_df is None or expansion_df.empty:
        return {}
    return _save_optional_csvs(directory, "Query genes", (
        ("query_genes_csv", f"{prefix}_query_genes.csv", expansion_df),
    ))


def save_proximity_results(directory: str, prefix: str, result: dict | None) -> dict:
    """把网络邻近度检验结果保存到分析目录({prefix}_proximity.csv,metric/value 两列)。报告页只作展示,这份文件是
    邻近度结果的唯一持久记录(cache 里的 JSON 按有效期自动清理);失败/未计算/空结果不落盘。约定同前。"""
    if not result or not result.get("ok"):
        return {}
    rows: list[dict] = []
    for key in ("d_c", "z", "p_empirical", "d_c_excl", "z_excl", "p_excl", "n_excl",
                "n_perm", "seed", "mean_random", "std_random", "lcc_nodes"):
        if result.get(key) is not None:
            rows.append({"metric": key, "value": result[key]})
    for key, val in (result.get("coverage") or {}).items():
        rows.append({"metric": f"coverage_{key}", "value": val})
    if not rows:
        return {}
    return _save_optional_csvs(directory, "Proximity results", (
        ("proximity_csv", f"{prefix}_proximity.csv", pd.DataFrame(rows)),
    ))
