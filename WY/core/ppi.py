"""PPI 网络分析:从 STRING 获取蛋白互作边表,构建网络并计算拓扑指标,绘图并保存结果。"""

import hashlib
import logging
import os

import pandas as pd
import requests

from .caching import open_csv_cache, save_csv_cache
from .http_client import request_with_retry
from .text_utils import safe_name


def _string_post(endpoint: str, params: dict) -> requests.Response | None:
    """向 STRING 数据库发送请求,失败时返回 None。
    提交的基因在 STRING 中一个都识别不出来(404),不算错误,按空结果处理,不重试。"""
    import config

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
    min_degree: int = 1,
) -> pd.DataFrame:
    """从 STRING 数据库获取基因列表的蛋白互作边表:先把基因符号映射为 STRING 蛋白 ID(不认识的基因被剔除),再按置信度阈值 score(0~1)筛选互作边。
    network_type 默认 "full" 为功能关联,"physical" 仅保留有物理结合证据的互作;species 为物种编号,人类是 "9606";基因列表自动去重排序,最多取 1000 个;min_degree 为建网后按连接数隐藏节点的阈值,默认 1 即只隐藏孤立节点(对齐文献惯例,不做更高连接度剪枝)。
    相同参数的结果存为本地缓存,有效期内直接复用;请求失败、基因全部不被识别或没有达到阈值的边时,返回空表。"""
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        return pd.DataFrame()
    if len(genes) > config.STRING_QUERY_LIMIT:
        genes = genes[:config.STRING_QUERY_LIMIT]

    # 缓存文件名由来源标识、置信度阈值、度数阈值拼成。
    # 相同参数组合固定同名,重新生成时直接覆盖旧文件。
    if label:
        label_part = safe_name(label)
        name = f"StringPPI_{label_part}_{score:g}_{int(min_degree)}.csv"
        if label_part == "manual":
            # 手动基因模式的基因列表无法从名字区分,附加基因列表短指纹作缓存标识,
            # 避免不同基因列表混用同一份缓存
            gene_fp = hashlib.md5("\r".join(genes).encode("utf-8")).hexdigest()[:8]
            name = f"StringPPI_manual_{score:g}_{int(min_degree)}_{gene_fp}.csv"
    else:
        name = f"StringPPI_{score:g}_{int(min_degree)}.csv"
    cache_file = os.path.join(config.PPI_CACHE_DIR, name)

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

    save_csv_cache(df, cache_file)
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


def render_core_target_plot(
    graph: "nx.Graph",
    top_genes: list[str],
    mcc_map: dict[str, float],
) -> "matplotlib.figure.Figure":
    """绘制 CytoHubba 风格的核心靶点子网络图(论文用):只画 top_genes 之间的互作,节点沿圆环按 MCC 排名
    从正上方顺时针排布,颜色从红(排名第 1)渐变到黄(排名末位),与 CytoHubba 默认配色一致,标签置于节点内部。"""
    import networkx as nx
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # 按排名顺序重建子图:circular_layout 按节点插入顺序布点,
    # 直接用 graph.subgraph 会保留原图顺序,排名与位置就对不上了
    sub = nx.Graph()
    sub.add_nodes_from(top_genes)
    for u, v, data in graph.edges(data=True):
        if u in sub and v in sub:
            sub.add_edge(u, v, **data)

    n = len(top_genes)
    pos = nx.circular_layout(sub)
    pos = {node: (-y, x) for node, (x, y) in pos.items()}  # 旋转 90°,第 1 名置于正上方

    cmap = cm.autumn  # 红(0) → 黄(1),与 CytoHubba 红黄渐变一致
    colors = [cmap(0.0)] if n == 1 else [cmap(i / (n - 1)) for i in range(n)]

    fig, ax = plt.subplots(figsize=(7.5, 7))
    nx.draw_networkx_edges(sub, pos, ax=ax, edge_color="#555555",
                           width=1.1, alpha=0.85)
    nx.draw_networkx_nodes(sub, pos, ax=ax, node_size=2300, node_color=colors,
                           edgecolors="#7a1f1f", linewidths=0.8)
    name_len = max(len(g) for g in top_genes)
    font_size = 9 if name_len <= 5 else (8 if name_len <= 7 else 7)
    nx.draw_networkx_labels(sub, pos, labels={g: g for g in top_genes},
                            font_size=font_size, font_weight="bold", ax=ax)
    ax.set_title(f"Core targets ranked by MCC (CytoHubba-style, top {n})",
                 fontsize=13)
    ax.axis("off")
    fig.tight_layout()
    return fig


# 模块着色用的 12 种颜色,循环使用
_MODULE_PALETTE = [
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948",
    "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
]


def render_pyvis_html(
    graph: "nx.Graph",
    metrics_df: pd.DataFrame,
    module_map: dict | None = None,
) -> str:
    """生成交互式 PPI 网络网页(HTML 文本):节点大小代表连接数多少,连线粗细代表置信度高低;默认 Hub 基因金色、其余节点淡蓝。
    传入 module_map(基因→模块ID)时改为按模块着色,不属于任何多基因模块的节点为浅灰;鼠标悬停节点可查看度数、模块等信息。"""
    import config
    from pyvis.network import Network

    hub_set = set(metrics_df.loc[metrics_df["is_hub"], "gene"])
    degree_map = dict(zip(metrics_df["gene"], metrics_df["degree"]))
    max_deg = max(metrics_df["degree"].max(), 1)

    module_color: dict[int, str] = {}
    if module_map:
        mids = sorted({int(v) for v in module_map.values() if int(v) > 0})
        module_color = {mid: _MODULE_PALETTE[(mid - 1) % len(_MODULE_PALETTE)] for mid in mids}

    nt = Network(height="600px", width="100%", directed=False, notebook=False)
    for node in graph.nodes():
        deg = degree_map.get(node, 0)
        size = 10 + int(40 * deg / max_deg)
        if module_map:
            mid = int(module_map.get(node, 0))
            color = module_color.get(mid, "#D9D9D9")  # 不在任何多基因模块中的节点
            title = f"{node}<br>degree: {deg}<br>module: {mid if mid > 0 else '-'}"
        else:
            color = "#F5C518" if node in hub_set else "#64B5F6"
            title = f"{node}<br>degree: {deg}"
        nt.add_node(
            node, label=node, size=size, color=color,
            title=title,
        )
    for u, v, data in graph.edges(data=True):
        s = float(data.get("score", 0))
        nt.add_edge(u, v, value=s / 1000, title=f"score: {s}")

    nt.set_options(f"""
    {{
      "physics": {{
        "solver": "forceAtlas2Based",
        "stabilization": {{"iterations": {config.PHYSICS_ITERATIONS}}},
        "forceAtlas2Based": {{
          "gravitationalConstant": -50,
          "springLength": 100,
          "springConstant": 0.08
        }}
      }},
      "interaction": {{"hover": true}}
    }}
    """)
    return nt.generate_html()


def render_static_plot(graph: "nx.Graph", metrics_df: pd.DataFrame) -> "matplotlib.figure.Figure":
    """用 networkx + matplotlib 绘制静态 PPI 网络图（论文用）:节点颜色映射 degree（附 colorbar）、大小 ∝ degree;只标注 Hub 基因标签避免拥挤。"""
    import networkx as nx
    import matplotlib.pyplot as plt
    from matplotlib import cm

    degrees = dict(graph.degree())
    node_degrees = [degrees[n] for n in graph.nodes()]
    # 交集过小且无互作边时所有节点度数全为 0,按 1 处理避免除零(节点统一用最小尺寸)
    max_deg = max(max(node_degrees, default=0), 1)
    n_nodes = max(len(graph.nodes()), 1)

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(graph, k=2.2 / (n_nodes ** 0.5), iterations=200, seed=42)

    nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.25, edge_color="#9e9e9e", width=0.8)
    nodes = nx.draw_networkx_nodes(
        graph, pos, ax=ax,
        node_size=[80 + 300 * d / max_deg for d in node_degrees],
        node_color=node_degrees,
        cmap=cm.viridis, vmin=0, vmax=max_deg,
        edgecolors="white", linewidths=0.6,
    )
    hub_set = set(metrics_df.loc[metrics_df["is_hub"], "gene"])
    nx.draw_networkx_labels(
        graph, pos, ax=ax,
        labels={n: n for n in graph.nodes() if n in hub_set}, font_size=8,
    )
    cb = fig.colorbar(nodes, ax=ax, shrink=0.8)
    cb.set_label("Degree")
    ax.set_title("PPI Network (node color = degree)", fontsize=14)
    ax.axis("off")
    fig.tight_layout()
    return fig


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
) -> dict:
    """把一次 PPI 分析的全部结果保存到结果目录下的 <prefix>_data 文件夹(prefix 为文件名前缀,如 aspirin_cancer),返回各文件的保存路径(含目录路径)。
    保存前会删除同 prefix 的旧文件夹(含历史结果),同一药物+疾病组合只保留最新一次分析,不会越积越多;节点表含 is_hub 列(算出 MCC 时含 mcc 列),其余为边表、交互式网络、GraphML、静态图矢量版及各富集/核心靶点/通路网络图表。
    传入富集结果表/图时才保存对应的富集文件,不传则跳过。"""
    import glob
    import shutil
    import networkx as nx
    import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    for old in glob.glob(os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data*")):
        if os.path.isdir(old):
            shutil.rmtree(old)
    subdir = os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data")
    os.makedirs(subdir, exist_ok=True)

    paths = {
        "nodes": os.path.join(subdir, f"{prefix}_nodes.csv"),
        "edges": os.path.join(subdir, f"{prefix}_edges.csv"),
        "html": os.path.join(subdir, f"{prefix}_network.html"),
        "graphml": os.path.join(subdir, f"{prefix}_graph.graphml"),
        "svg": os.path.join(subdir, f"{prefix}_network.svg"),
    }
    nodes_df.to_csv(paths["nodes"], index=False, encoding="utf-8-sig")
    edges_df.to_csv(paths["edges"], index=False, encoding="utf-8-sig")
    with open(paths["html"], "w", encoding="utf-8") as f:
        f.write(html_content)
    nx.write_graphml(graph, paths["graphml"])
    static_fig.savefig(paths["svg"], format="svg", bbox_inches="tight")

    if enrichment_df is not None and not enrichment_df.empty:
        paths["enrichment_csv"] = os.path.join(subdir, f"{prefix}_enrichment.csv")
        enrichment_df.to_csv(paths["enrichment_csv"], index=False, encoding="utf-8-sig")
    if enrichment_bar_fig is not None:
        paths["enrichment_bar_png"] = os.path.join(subdir, f"{prefix}_enrichment_bar.png")
        paths["enrichment_bar_svg"] = os.path.join(subdir, f"{prefix}_enrichment_bar.svg")
        enrichment_bar_fig.savefig(paths["enrichment_bar_png"], format="png", bbox_inches="tight")
        enrichment_bar_fig.savefig(paths["enrichment_bar_svg"], format="svg", bbox_inches="tight")
    if enrichment_dot_fig is not None:
        paths["enrichment_dot_png"] = os.path.join(subdir, f"{prefix}_enrichment_dot.png")
        paths["enrichment_dot_svg"] = os.path.join(subdir, f"{prefix}_enrichment_dot.svg")
        enrichment_dot_fig.savefig(paths["enrichment_dot_png"], format="png", bbox_inches="tight")
        enrichment_dot_fig.savefig(paths["enrichment_dot_svg"], format="svg", bbox_inches="tight")
    if core_fig is not None:
        paths["core_png"] = os.path.join(subdir, f"{prefix}_core_targets.png")
        paths["core_svg"] = os.path.join(subdir, f"{prefix}_core_targets.svg")
        core_fig.savefig(paths["core_png"], format="png", bbox_inches="tight")
        core_fig.savefig(paths["core_svg"], format="svg", bbox_inches="tight")
    if target_pathway_fig is not None:
        paths["dtp_png"] = os.path.join(subdir, f"{prefix}_enrichment_network.png")
        paths["dtp_svg"] = os.path.join(subdir, f"{prefix}_enrichment_network.svg")
        target_pathway_fig.savefig(paths["dtp_png"], format="png", bbox_inches="tight")
        target_pathway_fig.savefig(paths["dtp_svg"], format="svg", bbox_inches="tight")

    paths["directory"] = subdir
    return paths


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


def save_gsea_results(
    directory: str,
    prefix: str,
    drug_df: pd.DataFrame | None = None,
    disease_df: pd.DataFrame | None = None,
) -> dict:
    """把 GSEA prerank 结果追加保存到既有分析目录(<prefix>_data):包括药物靶点排秩结果({prefix}_gsea_drug.csv)
    和疾病靶点排秩结果({prefix}_gsea_disease.csv)。约定同前:空表不落盘。"""
    return _save_optional_csvs(directory, "GSEA results", (
        ("gsea_drug_csv", f"{prefix}_gsea_drug.csv", drug_df),
        ("gsea_disease_csv", f"{prefix}_gsea_disease.csv", disease_df),
    ))


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
