"""PPI 网络模块：STRING 边表抓取、networkx 建图与拓扑指标、PyVis/静态图渲染、结果落盘。"""

import hashlib
import logging
import os

import pandas as pd
import requests

from .caching import _is_cache_fresh, _read_cache_csv
from .text_utils import safe_name


def _string_post(endpoint: str, params: dict) -> requests.Response | None:
    """对 STRING API 发送 POST 请求，带超时与重试；最终失败返回 None。

    STRING 对"提交的标识符全部无法解析"的请求返回 404，这里视为空结果而非错误。
    """
    import time
    import config

    url = f"{config.STRING_BASE_URL}/{endpoint}"
    for attempt in range(config.API_MAX_RETRIES):
        try:
            resp = requests.post(url, data=params, timeout=config.API_TIMEOUT)
            if resp.status_code == 404:
                logging.warning("STRING %s: 404 (no identifiers resolved)", endpoint)
                return None
            resp.raise_for_status()
            return resp
        except Exception:
            logging.exception(
                "STRING %s request failed (attempt %d/%d)",
                endpoint, attempt + 1, config.API_MAX_RETRIES,
            )
            if attempt < config.API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)
    return None


def fetch_string_ppi(
    genes: list[str],
    species: str = "9606",
    score: float = 0.7,
    network_type: str = "full",
    label: str | None = None,
    min_degree: int = 3,
) -> pd.DataFrame:
    """
    通过 STRING API 获取基因列表的蛋白互作边表。

    请求逻辑分两步（STRING 官方最佳实践）：
      1. get_string_ids：把基因符号映射为 STRING 蛋白 ID，同时过滤掉
         STRING 不认识的基因（映射成功的才算数）；
      2. network：用映射后的 STRING ID 获取相互作用边，按 required_score
         （= score×1000）过滤低置信度连接。

    结果缓存在 Static/data/cache/ppi 下，文件名为
    StringPPI_{label}_{STRING confidence}_{min_degree}.csv
    （label 为来源标识，如 药物名_疾病名；无 label 时退化为 StringPPI_{confidence}_{min_degree}.csv）。
    手动基因模式（label="manual"）额外附加基因列表短指纹（md5 前 8 位）参与键。
    相同参数组合固定同名，重新生成时直接覆盖旧文件，不在缓存目录堆积；
    TTL 内命中缓存直接读本地，不访问远端 API。

    Parameters
    ----------
    genes : list[str]
        基因符号列表（如 ["TP53", "BRCA1"]），会自动去重排序并截断到 1000。
    species : str
        NCBI taxon id，如 "9606"（人类）。
    score : float
        置信度阈值（0~1），映射为 STRING 的 required_score（0~1000）。
    network_type : str
        "full"（功能关联，默认）或 "physical"（仅物理结合）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。
    min_degree : int
        Hub 筛选的度数阈值（mark_hub_genes 用），仅参与缓存文件命名，默认 3。

    Returns
    -------
    pd.DataFrame
        每行一条相互作用边，含 protein_A / protein_B / score 及各证据分数列；
        请求失败、全部基因不被识别或没有符合阈值的边时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        return pd.DataFrame()
    if len(genes) > config.STRING_QUERY_LIMIT:
        genes = genes[:config.STRING_QUERY_LIMIT]

    os.makedirs(config.PPI_CACHE_DIR, exist_ok=True)
    # 命名规则：StringPPI_{药物名_疾病名}_{STRING confidence}_{Min degree}.csv
    # 相同参数组合固定同名，重新生成时 to_csv 直接覆盖旧文件。
    if label:
        label_part = safe_name(label)
        name = f"StringPPI_{label_part}_{score:g}_{int(min_degree)}.csv"
        if label_part == "manual":
            # 手动基因模式：基因列表不由药物/疾病名决定，附加短指纹（md5 前 8 位）
            # 参与缓存键，防止 TTL 内不同基因列表互相污染缓存
            gene_fp = hashlib.md5("\r".join(genes).encode("utf-8")).hexdigest()[:8]
            name = f"StringPPI_manual_{score:g}_{int(min_degree)}_{gene_fp}.csv"
    else:
        name = f"StringPPI_{score:g}_{int(min_degree)}.csv"
    cache_file = os.path.join(config.PPI_CACHE_DIR, name)

    # ---- cache hit: load local (若超过 TTL 视为过期，重新请求) ----
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached STRING PPI for %d genes", len(genes))
        return _read_cache_csv(cache_file)

    # ---- Step 1: 标识符映射（基因符号 → STRING ID） ----
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

    # ---- Step 2: 获取相互作用边 ----
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

    # 只保留输入基因内部的边（大小写不敏感），整理列名
    gene_set = {g.upper() for g in genes}
    df = df[
        df["preferredName_A"].str.upper().isin(gene_set)
        & df["preferredName_B"].str.upper().isin(gene_set)
    ].rename(columns={"preferredName_A": "protein_A", "preferredName_B": "protein_B"})

    keep = ["protein_A", "protein_B", "score", "nscore", "fscore", "pscore",
            "ascore", "escore", "dscore", "tscore"]
    df = df[[c for c in keep if c in df.columns]]

    # ---- cache to CSV ----
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("STRING PPI for %d genes: %d edges cached", len(genes), len(df))
    return df


def build_ppi_graph(edges_df: pd.DataFrame, all_genes: list[str]) -> "nx.Graph":
    """
    用 STRING 边表构建 networkx 无向图。

    - 全部输入基因都作为节点（无边的基因为孤立节点，degree=0）；
    - 去除自环与重复边（同一对基因保留置信度最高的一条）；
    - 边附带 score（STRING 原始分 0~1000）与 weight（归一化 0~1）。
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
    计算网络节点的拓扑指标：degree / betweenness / closeness。

    孤立节点（degree=0）的中心性为 0；结果按 degree 降序排列。
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
    标记 Hub 基因：degree 大于全部节点 degree 的中位数，且 degree ≥ min_degree（兜底阈值）。

    返回原表并新增 is_hub 布尔列。
    """
    df = metrics_df.copy()
    threshold = df["degree"].median()
    df["is_hub"] = (df["degree"] > threshold) & (df["degree"] >= min_degree)
    return df


def render_pyvis_html(graph: "nx.Graph", metrics_df: pd.DataFrame) -> str:
    """
    用 PyVis 生成交互式 PPI 网络 HTML。

    - 节点大小 ∝ degree，Hub 基因金色高亮、普通节点淡蓝；
    - 边宽 ∝ 置信度分数；
    - 物理引擎迭代次数取 config.PHYSICS_ITERATIONS，防止节点过度拥挤。

    Returns
    -------
    str
        HTML 内容字符串（可直接用 st.iframe 内嵌或保存为 .html 文件）。
    """
    import config
    from pyvis.network import Network

    hub_set = set(metrics_df.loc[metrics_df["is_hub"], "gene"])
    degree_map = dict(zip(metrics_df["gene"], metrics_df["degree"]))
    max_deg = max(metrics_df["degree"].max(), 1)

    nt = Network(height="600px", width="100%", directed=False, notebook=False)
    for node in graph.nodes():
        deg = degree_map.get(node, 0)
        size = 10 + int(40 * deg / max_deg)
        color = "#F5C518" if node in hub_set else "#64B5F6"
        nt.add_node(
            node, label=node, size=size, color=color,
            title=f"{node}<br>degree: {deg}",
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
    将一次 PPI 分析的全部结果保存到 DOWNLOADS_DIR/<prefix>_data/ 目录（替代式）。

    目录名 = <prefix>_data（prefix 如 药物名_疾病名，例如 aspirin_cancer_data），
    文件名 = <prefix>_xxx（如 aspirin_cancer_nodes.csv）。

    【替代式保存】：保存前会删除同 prefix 的旧目录（含历史时间戳目录），
    同一 prefix（同一药物+疾病）始终只保留最新一次分析，不会堆积。

    目录内容：
      - {prefix}_nodes.csv      节点拓扑指标表（含 is_hub 列）
      - {prefix}_edges.csv      边表（protein_A / protein_B / score ...）
      - {prefix}_network.html   交互式网络
      - {prefix}_graph.graphml  可编辑网络文件（供 Cytoscape / Gephi 精修）
      - {prefix}_network.svg    静态图矢量版（论文出版用）
      - {prefix}_enrichment.csv           富集结果表（有富集数据时）
      - {prefix}_enrichment_bar.png/.svg  富集柱状图（有富集数据时）
      - {prefix}_enrichment_dot.png/.svg  富集气泡图（有富集数据时）

    Parameters
    ----------
    graph : networkx.Graph
        分析构建的网络图，用于 GraphML 导出。
    static_fig : matplotlib.figure.Figure
        静态网络图，用于 SVG 导出。
    nodes_df : pd.DataFrame
        节点拓扑指标表。
    edges_df : pd.DataFrame
        边表。
    html_content : str
        交互式网络的 HTML 内容。
    prefix : str
        文件名/目录名前缀，如 "aspirin_cancer"（目录会追加 "_data" 后缀）。
    enrichment_df : pd.DataFrame | None
        富集结果表（统一列）；None 表示本次分析无富集数据。
    enrichment_bar_fig / enrichment_dot_fig : Figure | None
        富集柱状图 / 气泡图；None 表示无富集数据。

    Returns
    -------
    dict
        各文件路径与目录路径（含可选富集文件路径）。
    """
    import glob
    import shutil
    import networkx as nx
    import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    # 替代式保存：删除同前缀的旧目录（含历史时间戳目录），只保留本次最新分析
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

    # 富集文件（可选）
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
