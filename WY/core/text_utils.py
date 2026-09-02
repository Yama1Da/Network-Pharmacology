"""通用文本 / 数值小工具。"""


def safe_name(name: str) -> str:
    """把用户输入(药物名/疾病名/基因等)转成能安全用作文件名的形式。

    把文件名不允许出现的符号统一替换成下划线,防止输入内容被当作路径、
    破坏文件组织;转换结果为空时用 "unnamed"。只用于文件命名,
    查询时仍使用原始输入。
    """
    import re

    cleaned = re.sub(r"[^\w\-]+", "_", str(name).strip())
    return cleaned.strip("._-") or "unnamed"


def _safe_round(value, ndigits: int = 4) -> float:
    """安全取数:部分靶点在数据库中没有指标数值,统一按 0 处理后再四舍五入。"""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return 0.0
