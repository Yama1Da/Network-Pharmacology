"""PPI 网络图渲染器。

核心靶点环图(CytoHubba 风格)、pyvis 交互网页、静态网络图,
以及 STRING / Cytoscape 两种论文风格图。全部为纯渲染函数:
图数据(networkx 图、指标表、边表)由调用方传入,不取数、不落盘。
"""

import functools

import pandas as pd


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
    from WY import config
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


# STRING 官方证据通道配色:列名 → (颜色, 图例名);连线颜色取七个通道中得分最高者
_EVIDENCE_CHANNELS = (
    ("escore", "#F963C1", "Experimentally determined"),
    ("dscore", "#41CFC0", "Curated databases"),
    ("tscore", "#E8E33B", "Textmining"),
    ("ascore", "#111111", "Coexpression"),
    ("pscore", "#3E5BD9", "Co-occurrence"),
    ("nscore", "#2FBF3F", "Gene neighborhood"),
    ("fscore", "#E63B3B", "Gene fusion"),
)


def _edge_evidence_map(edges_df: pd.DataFrame) -> dict:
    """给每条互作边标注主导证据通道:取七个 STRING 通道得分里最高的那个。
    返回 {(u, v): (颜色, 图例名, 得分)}(键为蛋白对排序后的元组);通道列缺失或得分全零时返回空表,
    调用方据此把所有边画成灰色。"""
    # 实际存在的通道列单独建档:边表可能缺个别通道列(旧缓存/API 字段调整),
    # 主导通道必须按"过滤后"的列表回查,直接下标 _EVIDENCE_CHANNELS 会整体错位
    channel_meta = [(c, color, legend) for c, color, legend in _EVIDENCE_CHANNELS
                    if c in edges_df.columns]
    if not channel_meta or edges_df.empty:
        return {}
    channel_cols = [c for c, _, _ in channel_meta]
    vals = edges_df[channel_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    for c in channel_cols:  # 兼容个别通道以 0-999 计数返回的情况
        if float(vals[c].max()) > 1.5:
            vals[c] = vals[c] / 1000.0

    mapping: dict[tuple[str, str], tuple[str, str, float]] = {}
    order = edges_df.sort_values("score", ascending=False).index
    seen: set[tuple[str, str]] = set()
    for idx in order:
        row = vals.loc[idx]
        u = str(edges_df.at[idx, "protein_A"])
        v = str(edges_df.at[idx, "protein_B"])
        key = tuple(sorted((u, v)))
        if key in seen:
            continue
        seen.add(key)
        arr = row.to_numpy(dtype=float)
        best = int(arr.argmax())
        if arr[best] <= 0:
            continue  # 七个通道全零:没有主导证据,不进映射,调用方按灰色画
        _, color, legend = channel_meta[best]
        mapping[key] = (color, legend, float(arr[best]))
    return mapping


@functools.lru_cache(maxsize=8)
def _ball_sprite(color: str, px: int = 72) -> tuple:
    """预渲染一张径向渐变"蛋白球"贴图(基色 + 左上高光 + 边缘压暗),模拟 STRING 的 3D 球节点。
    返回 (rgb 数组, alpha 数组),均为 px×px;按颜色缓存,同色节点复用同一张贴图。"""
    import numpy as np
    from matplotlib import colors as mcolors

    yy, xx = np.mgrid[0:px, 0:px]
    u = xx / (px - 1) * 2.0 - 1.0
    v = yy / (px - 1) * 2.0 - 1.0
    r = np.hypot(u, v)
    base = np.array(mcolors.to_rgb(color))
    shade = 0.38 + 0.62 * np.clip(1.0 - r, 0.0, 1.0) ** 0.85
    spec = np.exp(-(((u + 0.42) ** 2 + (v + 0.42) ** 2) / 0.14))
    rgb = np.clip(base * shade[..., None] + spec[..., None] * 0.85, 0.0, 1.0)
    alpha = np.clip((1.0 - r) * px * 0.5, 0.0, 1.0)
    return rgb, alpha


def render_string_style_plot(
    graph: "nx.Graph",
    metrics_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> "matplotlib.figure.Figure":
    """绘制 STRING 风格的静态 PPI 网络图(论文用):节点为带高光的渐变"蛋白球",大小 ∝ 连接数,
    Hub 基因金色、其余淡蓝,标签置于节点外侧;连线颜色 = 该互作最强的证据通道(STRING 官方配色,
    左下角附图例),粗细 ∝ 置信度。布局沿用固定 seed 的 spring 布局保证可复现。
    edges_df 为 STRING 原始边表,用于提取各证据通道得分;通道列缺失或全零时连线退为灰色。"""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import networkx as nx
    import numpy as np

    n_nodes = len(graph.nodes())

    fig, ax = plt.subplots(figsize=(13, 11))
    if n_nodes == 0:
        ax.text(0.5, 0.5, "No interactions above threshold", ha="center",
                va="center", fontsize=12, color="#777777")
        ax.axis("off")
        return fig

    hub_set = (set(metrics_df.loc[metrics_df["is_hub"], "gene"])
               if "is_hub" in metrics_df.columns else set())
    degree_map = dict(zip(metrics_df["gene"], metrics_df["degree"]))
    max_deg = max(int(metrics_df["degree"].max()), 1) if len(metrics_df) else 1

    pos = nx.spring_layout(graph, k=2.2 / (n_nodes ** 0.5), iterations=200, seed=42)
    radius = {n: 0.014 + 0.030 * (float(degree_map.get(n, 0) / max_deg) ** 0.5)
              for n in graph.nodes()}

    # 边:先画灰底衬托,再按主导证据通道着色(通道缺失/全零时保持灰色)
    evidence = _edge_evidence_map(edges_df)
    edge_colors: list[str] = []
    edge_widths: list[float] = []
    for u, v, data in graph.edges(data=True):
        color, _, best = evidence.get(tuple(sorted((str(u), str(v)))),
                                      ("#9a9a9a", "", float(data.get("score", 0)) / 1000))
        edge_colors.append(color)
        edge_widths.append(0.4 + 1.8 * min(max(best, 0.0), 1.0))
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color=edge_colors,
                           width=edge_widths, alpha=0.8)

    # 节点:渐变球贴图按圆形裁剪逐个贴上
    for node, (x, y) in pos.items():
        r = radius[node]
        color = "#F5C518" if node in hub_set else "#64B5F6"
        rgb, alpha = _ball_sprite(color)
        im = ax.imshow(rgb, extent=(x - r, x + r, y - r, y + r), zorder=3,
                       interpolation="bilinear", origin="upper")
        im.set_alpha(alpha)
        im.set_clip_path(mpatches.Circle((x, y), r, transform=ax.transData))

    # 标签:沿"布局重心 → 节点"的离心方向外推到节点外,字号随网络规模收缩
    cx = float(np.mean([p[0] for p in pos.values()]))
    cy = float(np.mean([p[1] for p in pos.values()]))
    font_by_scale = 7 if n_nodes <= 60 else (6 if n_nodes <= 120 else
                                             (5 if n_nodes <= 200 else 4))
    for node, (x, y) in pos.items():
        dx, dy = x - cx, y - cy
        dist = float(np.hypot(dx, dy)) or 1.0
        lx = x + dx / dist * (radius[node] + 0.02)
        ly = y + dy / dist * (radius[node] + 0.02)
        ax.text(lx, ly, node, fontsize=font_by_scale,
                fontweight="bold" if node in hub_set else "normal",
                ha="center", va="center", color="#1a1a1a", zorder=4)

    # 显式锁定坐标范围:imshow 会触发自动缩放,把窗口塌缩到最后一张贴图的范围内,
    # 其余元素被甩出画布(tight bbox 随之爆炸),必须在此收口
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.15, max(xs) + 0.15)
    ax.set_ylim(min(ys) - 0.15, max(ys) + 0.15)
    ax.set_aspect("equal", adjustable="box")

    # 图例:只列本次网络里实际出现的证据通道
    present = {legend for _, legend, _ in evidence.values()}
    if present:
        handles = [Line2D([0], [0], color=color, lw=2.8)
                   for _, color, legend in _EVIDENCE_CHANNELS if legend in present]
        labels = [legend for _, _, legend in _EVIDENCE_CHANNELS if legend in present]
        ax.legend(handles, labels, loc="lower left", fontsize=7,
                  title="Interaction evidence", title_fontsize=8,
                  framealpha=0.9, borderpad=0.6)

    ax.set_title("STRING-style PPI network (edge color = interaction evidence)",
                 fontsize=14)
    ax.axis("off")
    return fig


def render_cytoscape_style_plot(graph: "nx.Graph", metrics_df: pd.DataFrame) -> "matplotlib.figure.Figure":
    """绘制 Cytoscape 风格(论文常用)的静态 PPI 网络图:节点按连接数从高到低依次排入同心圆环,
    高连接数的 Hub 居中、低连接数靶点在外圈;节点大小 ∝ 连接数,颜色按 viridis 反向渐变
    (高 = 深紫,低 = 黄),标签置于节点内部并随节点尺寸自适应;连线为细灰边。布局完全由度数
    决定(确定性,无随机力导向),同环节点按周长约束外推环半径避免重叠。"""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    import networkx as nx
    import numpy as np

    deg_map = dict(zip(metrics_df["gene"], metrics_df["degree"].astype(int)))
    nodes = sorted(graph.nodes(), key=lambda v: (-int(deg_map.get(v, 0)), str(v)))
    n = len(nodes)

    fig, ax = plt.subplots(figsize=(13, 12))
    if n == 0:
        ax.text(0.5, 0.5, "No nodes to display", ha="center", va="center",
                fontsize=12, color="#777777")
        ax.axis("off")
        return fig

    dmax = max(int(deg_map.get(nodes[0], 0)), 1)
    radius = {v: 0.10 + 0.24 * float(int(deg_map.get(v, 0)) / dmax) ** 0.5
              for v in nodes}

    # 环容量按周长比例分配(内圈小、外圈大),舍入误差全部并入最外圈
    n_rings = max(2, min(8, int(np.ceil((n ** 0.5) / 2.2))))
    weights = np.arange(1, n_rings + 1, dtype=float)
    caps = np.maximum(1, np.round(n * weights / weights.sum()).astype(int))
    diff = int(caps.sum()) - n
    i = n_rings - 1
    while diff > 0:  # 舍入超出:从外圈往内逐环回收名额
        if caps[i] > 1:
            caps[i] -= 1
            diff -= 1
        i = n_rings - 1 if i == 0 else i - 1
    if diff < 0:
        caps[n_rings - 1] += -diff

    rings: list[list[str]] = []
    start = 0
    for cap in caps:
        rings.append(nodes[start:start + int(cap)])
        start += int(cap)
    rings = [r for r in rings if r]

    # 环半径:相邻环留出"双方最大节点半径 + 间隙",且周长须容纳本环全部节点
    ring_R: list[float] = []
    prev_R = prev_r = 0.0
    for i, ring in enumerate(rings):
        r_i = max(radius[v] for v in ring)
        if i == 0 and len(ring) == 1:
            R = 0.0  # 唯一的第 1 名放在圆心(与参考图一致)
        else:
            R = prev_R + prev_r + r_i + 0.30
            R = max(R, (2 * r_i + 0.30) * len(ring) / (2 * np.pi))
        ring_R.append(R)
        prev_R, prev_r = R, r_i

    pos: dict[str, tuple[float, float]] = {}
    for i, ring in enumerate(rings):
        R, k = ring_R[i], len(ring)
        if R == 0.0:
            pos[ring[0]] = (0.0, 0.0)
            continue
        for j, node in enumerate(ring):
            ang = np.pi / 2 + (j + (0.5 if i % 2 else 0.0)) * 2 * np.pi / k
            pos[node] = (R * float(np.cos(ang)), R * float(np.sin(ang)))

    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#9a9a9a",
                           width=0.7, alpha=0.45)

    cmap_r = plt.get_cmap("viridis_r")
    norm = Normalize(vmin=0, vmax=dmax)
    for node, (x, y) in pos.items():
        ax.add_patch(mpatches.Circle(
            (x, y), radius[node],
            facecolor=cmap_r(norm(int(deg_map.get(node, 0)))),
            edgecolor="#3A3A3A", linewidth=0.5, zorder=3))

    # 标签置于节点内部:字号按节点半径与基因名长度自适应(数据单位→pt 换算后钳制)
    outer_R = ring_R[-1] + max(radius[v] for v in rings[-1])
    ax.set_xlim(-outer_R - 1.0, outer_R + 1.0)
    ax.set_ylim(-outer_R - 1.0, outer_R + 1.0)
    ax.set_aspect("equal")
    pt_per_unit = 72.0 * 13.0 / (2 * (outer_R + 1.0))
    for node, (x, y) in pos.items():
        fs = float(np.clip(2 * radius[node] * pt_per_unit * 0.80
                           / (0.60 * max(len(node), 3)), 3.2, 9.5))
        # 深色节点(高连接数,近紫色)配白字,浅色节点配黑字,保证标签可读
        fc = cmap_r(norm(int(deg_map.get(node, 0))))[:3]
        lum = 0.299 * fc[0] + 0.587 * fc[1] + 0.114 * fc[2]
        ax.text(x, y, node, fontsize=fs, ha="center", va="center",
                color="#ffffff" if lum < 0.45 else "#111111", zorder=4)

    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap_r),
                      ax=ax, shrink=0.55, pad=0.02)
    cb.set_label("Degree", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    ax.set_title("PPI network (Cytoscape-style, degree-sorted concentric rings)",
                 fontsize=14)
    ax.axis("off")
    return fig
