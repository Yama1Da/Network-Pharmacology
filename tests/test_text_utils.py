"""core.text_utils 的测试:文件名消毒、安全取数、基因列表解析、名称映射回环。

阶段三第 1 批新写(原自检脚本未覆盖此模块)。
"""
from WY.core.text_utils import (
    genes_from_dataframe,
    load_name_map,
    parse_genes_text,
    remember_name_mapping,
    resolve_original_name,
    safe_name,
    safe_round,
)


def test_safe_name_symbols_and_edges():
    """文件名消毒:符号归下划线、空串兜底 unnamed、路径分隔符不外逃。"""
    assert safe_name("metformin") == "metformin"          # 无损名原样
    assert safe_name("  ") == "unnamed"                    # 空串兜底
    assert safe_name("a/b\\c:d*e") == "a_b_c_d_e"          # 路径分隔符与符号
    assert safe_name("St. John's Wort") == "St_John_s_Wort"  # 含点/撇号(与对接归档前缀一致)
    assert safe_name("二甲双胍") == "二甲双胍"               # 中文保留


def test_safe_name_windows_reserved():
    """Windows 保留设备名(CON/NUL/COM1…)加下划线前缀,否则建文件必失败。"""
    assert safe_name("CON") == "_CON"
    assert safe_name("nul") == "_nul"
    assert safe_name("COM1") == "_COM1"
    assert safe_name("CONNECTION") == "CONNECTION"  # 只是前缀撞上,不算保留名


def test_safe_round():
    """安全取数:正常四舍五入;None/非数值按 0.0;数字字符串可转换。"""
    assert safe_round(3.14159, 4) == 3.1416
    assert safe_round("3.14159", 2) == 3.14
    assert safe_round(None) == 0.0
    assert safe_round("abc") == 0.0


def test_parse_genes_text_formats():
    """基因列表字段解析:列表样式/分号/逗号/空白全兼容;空输入返回空列表。"""
    assert parse_genes_text("['TP53', 'AKT1', 'MYC']") == ["TP53", "AKT1", "MYC"]
    assert parse_genes_text("TP53;STAT3;TP53") == ["TP53", "STAT3", "TP53"]
    assert parse_genes_text("IL6,PTGS2") == ["IL6", "PTGS2"]
    assert parse_genes_text("TP53 AKT1") == ["TP53", "AKT1"]
    assert parse_genes_text("") == []
    assert parse_genes_text(None) == []


def test_genes_from_dataframe():
    """汇总富集表各条目基因:统一大写、跨行去重。"""
    import pandas as pd

    df = pd.DataFrame({
        "genes": ["['tp53', 'akt1']", "TP53;VEGFA", ""],
    })
    assert genes_from_dataframe(df) == {"TP53", "AKT1", "VEGFA"}


def test_name_mapping_roundtrip():
    """名称映射:有损转换记映射并可还原;无损转换不写盘;缺记录返回 None。"""
    remember_name_mapping("Type 2 Diabetes")     # 有损(空格→下划线)→ 记录
    key = safe_name("Type 2 Diabetes")
    assert load_name_map().get(key) == "Type 2 Diabetes"
    assert resolve_original_name(key) == "Type 2 Diabetes"
    # 无损名字不写映射
    remember_name_mapping("metformin")
    assert "metformin" not in load_name_map()
    assert resolve_original_name("metformin") is None
    # 缺记录返回 None
    assert resolve_original_name("NOPE_NAME_X") is None
