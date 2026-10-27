"""core.geo 的测试:样本分组(丢弃剔除与四道硬边界)、样本来源启发式标注。

迁移自 selfcheck_ml.py(17 项中的 3 项)。
"""
import pandas as pd

from WY.core.geo import GeoDataset, assign_groups, sample_source_label


def _synth_samples(groups: tuple) -> GeoDataset:
    """按 (组名, 样本数) 列表造一个只有样本注释、无表达矩阵的合成 GeoDataset。"""
    rows = []
    for g, n in groups:
        for i in range(n):
            rows.append((f"GSM_{g}_{i}", "tissue sample", f"group: {g}"))
    samples = pd.DataFrame(
        {"source_name_ch1": [r[1] for r in rows],
         "characteristics_ch1": [r[2] for r in rows]},
        index=[r[0] for r in rows])
    return GeoDataset(gse="SYNTH", meta={}, samples=samples,
                      table=pd.DataFrame(), gpl="")


def test_dropped_samples():
    """M1 审查回归:未命中任何 pattern 的样本应从标签中真正剔除。"""
    ds = _synth_samples((("Severe", 5), ("Incipient", 4), ("Moderate", 4), ("Control", 6)))
    labels = assign_groups(ds, "Severe", "Control")
    n_case, n_ctrl, n_total = int((labels == 1).sum()), int((labels == 0).sum()), len(labels)
    assert n_case == 5 and n_ctrl == 6, f"分组数错: case={n_case} ctrl={n_ctrl}"
    assert n_total == 11, f"丢弃样本未剔除: 总数 {n_total}(期望 11 = 5+6)"
    assert "GSM_Incipient_0" not in labels.index, "Incipient 样本不应在标签中"
    assert "GSM_Moderate_0" not in labels.index, "Moderate 样本不应在标签中"


def test_group_hard_boundaries():
    """分组四道硬边界:零命中/全命中/重叠/样本不足均应报错。"""
    ds = _synth_samples((("Severe", 5), ("Incipient", 4), ("Control", 6)))

    for case_p, ctrl_p, expect in [
        ("zzz", "Control", "零命中"),
        ("tissue", "Control", "全命中"),
        ("Severe|Control", "Control", "重叠"),
        ("zzz", "Control", "零命中(case 匹配 0 个)"),
    ]:
        try:
            assign_groups(ds, case_p, ctrl_p)
            raise AssertionError(f"{expect}: 未报错")
        except RuntimeError:
            pass  # 预期报错


def test_sample_source_label():
    """样本来源启发式(方案 §4.1/§4.3):检索表 possible_cellline 与 shortlist
    evidence_scope 共用 sample_source_label,三值口径 + 细胞系优先。"""
    assert sample_source_label("HEPG2 cell line", "hepatic cells in vitro") == "细胞系", \
        "细胞系关键词未命中"
    assert sample_source_label("liver tissue biopsy", "patients") == "组织", \
        "组织关键词未命中"
    assert sample_source_label("", "a randomized double-blind trial") == "未标注", \
        "无关键词应判未标注"
    assert sample_source_label("cell line and tissue samples", "") == "细胞系", \
        "同时命中时细胞系应优先"
