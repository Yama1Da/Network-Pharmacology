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

    如果提交的基因在 STRING 中一个都识别不出来(404),不算错误,按空结果处理,不重试。
    """
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
    min_degree: int = 3,
) -> pd.DataFrame:
    """
    从 STRING 数据库获取基因列表的蛋白互作边表。

    分两步:先把基因符号映射为 STRING 蛋白 ID(STRING 不认识的基因被剔除),
    再按置信度阈值 score(0~1)筛选互作边。network_type 默认 "full" 为功能关联,
    "physical" 仅保留有物理结合证据的互作;species 为物种编号,人类是 "9606"。
    基因列表会自动去重排序,最多取 1000 个。

    相同参数的结果存为本地缓存,有效期内直接复用,不再重复请求。
    请求失败、基因全部不被识别或没有达到阈值的边时,返回空表。
    """
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

    # ---- 命中有效期内的缓存:直接读本地文件 ----
    cached = open_csv_cache(cache_file, fresh=True, label=f"STRING PPI for {len(genes)} genes")
    if cached is not None:
        return cached

    # ---- 第一步:把基因符号映射为 STRING 蛋白 ID ----
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

    # ---- 第二步:获取蛋白相互作用边 ----
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

    # 只保留两端都是输入基因的边,并统一列名
    gene_set = {g.upper() for g in genes}
    df = df[
        df["preferredName_A"].str.upper().isin(gene_set)
        & df["preferredName_B"].str.upper().isin(gene_set)
    ].rename(columns={"preferredName_A": "protein_A", "preferredName_B": "protein_B"})

    keep = ["protein_A", "protein_B", "score", "nscore", "fscore", "pscore",
            "ascore", "escore", "dscore", "tscore"]
    df = df[[c for c in keep if c in df.columns]]

    # ---- 结果存入本地缓存 ----
    save_csv_cache(df, cache_file)
    logging.info("STRING PPI for %d genes: %d edges cached", len(genes), len(df))
    return df


def build_ppi_graph(edges_df: pd.DataFrame, all_genes: list[str]) -> "nx.Graph":
    """
    用 STRING 边表构建蛋白互作网络图。

    所有输入基因都作为节点(没有互作关系的基因成为孤立节点);同一对基因的
    重复边只保留置信度最高的一条。
    """
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(all_genes)
    if edges_df.empty:
        return G

    edges = edges_df[edges_df["protein_A"] != edges_df["protein_B"]].copy()
    # 统一 A/B 方向后按键去重，保留 score 最高的边
    edges["_key"] = edges.apply(
        lambda r: tuple(sorted((str(r["protein_A"]), str(r["protein_B"])))), axis=1
    )
    edges = edges.sort_values("score", ascending=False).drop_duplicates(subset="_key")
    for _, r in edges.iterrows():
        s = float(r["score"])
        G.add_edge(r["protein_A"], r["protein_B"], score=s, weight=s / 1000.0)
    return G


def compute_ppi_metrics(graph: "nx.Graph") -> pd.DataFrame:
    """
    计算每个节点的拓扑指标:度(degree)、介数中心性、接近中心性。

    孤立节点各项为 0;结果按度数从高到低排列。
    """
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


def mark_hub_genes(metrics_df: pd.DataFrame, min_degree: int = 3) -> pd.DataFrame:
    """
    标记 Hub 基因:连接数高于全部节点中位数、且不低于 min_degree 的节点。

    在原表基础上新增 is_hub 一列(True 为 Hub 基因)。
    """
    df = metrics_df.copy()
    threshold = df["degree"].median()
    df["is_hub"] = (df["degree"] > threshold) & (df["degree"] >= min_degree)
    return df


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
    """
    生成交互式 PPI 网络网页(HTML 文本)。

    节点大小代表连接数多少,连线粗细代表置信度高低;默认 Hub 基因金色、
    其余节点淡蓝。传入 module_map(基因→模块ID)时改为按模块着色,不属于
    任何多基因模块的节点为浅灰;鼠标悬停节点可查看度数、模块等信息。
    """
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
    """
    用 networkx + matplotlib 绘制静态 PPI 网络图（论文用）。

    节点颜色映射 degree（附 colorbar）、大小 ∝ degree；只标注 Hub 基因标签避免拥挤。
    """
    import networkx as nx
    import matplotlib.pyplot as plt
    from matplotlib import cm

    degrees = dict(graph.degree())
    node_degrees = [degrees[n] for n in graph.nodes()]
    max_deg = max(node_degrees, default=1)
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
) -> dict:
    """
    把一次 PPI 分析的全部结果保存到结果目录下的 <prefix>_data 文件夹。

    保存前会删除同 prefix 的旧文件夹(含历史结果),同一药物+疾病组合只保留
    最新一次分析,不会越积越多。prefix 为文件名前缀(如 aspirin_cancer)。

    目录内容：
      - {prefix}_nodes.csv      节点拓扑指标表（含 is_hub 列）
      - {prefix}_edges.csv      边表（protein_A / protein_B / score ...）
      - {prefix}_network.html   交互式网络
      - {prefix}_graph.graphml  可编辑网络文件（供 Cytoscape / Gephi 精修）
      - {prefix}_network.svg    静态图矢量版（论文出版用）
      - {prefix}_enrichment.csv           富集结果表（有富集数据时）
      - {prefix}_enrichment_bar.png/.svg  富集柱状图（有富集数据时）
      - {prefix}_enrichment_dot.png/.svg  富集气泡图（有富集数据时）

    传入富集结果表/图时才保存对应的富集文件,不传则跳过。
    返回各文件的保存路径(含目录路径)。
    """
    import glob
    import shutil
    import networkx as nx
    import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    # 删除同前缀的旧结果文件夹,只保留本次分析
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

    # 富集结果文件(没有则不保存)
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

    paths["directory"] = subdir
    return paths


def _save_optional_csvs(directory: str, label: str,
                        specs: tuple[tuple[str, str, pd.DataFrame | None], ...]) -> dict:
    """把若干结果表追加保存到既有分析目录(<prefix>_data)。

    specs 每项为 (结果键, 文件名, 表格)；约定:空表不落盘,避免产生无意义的空文件。
    返回实际写出的文件路径,调用方会并入结果字典,打包下载时自动包含。
    """
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
    """
    把蛋白模块检测结果追加保存到既有分析目录(<prefix>_data)。

    包括模块汇总表({prefix}_modules.csv)和模块富集长表
    ({prefix}_module_enrichment.csv)。约定与 _save_optional_csvs 相同:空表不落盘。
    """
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
    """
    把 GSEA prerank 结果追加保存到既有分析目录(<prefix>_data)。

    包括药物靶点排秩结果({prefix}_gsea_drug.csv)和疾病靶点排秩结果
    ({prefix}_gsea_disease.csv)。约定同前:空表不落盘。
    """
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
    """
    把 RWR 靶点优先级排序结果追加保存到既有分析目录(<prefix>_data)。

    包括以疾病靶点为种子的排序表({prefix}_rwr_disease_seeds.csv)和以药物
    靶点为种子的排序表({prefix}_rwr_drug_seeds.csv)。约定同前:空表不落盘。
    """
    return _save_optional_csvs(directory, "RWR results", (
        ("rwr_disease_csv", f"{prefix}_rwr_disease_seeds.csv", disease_df),
        ("rwr_drug_csv", f"{prefix}_rwr_drug_seeds.csv", drug_df),
    ))
