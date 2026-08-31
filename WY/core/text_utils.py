"""通用文本 / 数值小工具。"""


def safe_name(name: str) -> str:
    """把用户输入（药物名/疾病名/基因等）净化为可安全用作文件名的字符串。

    将所有非字母/数字/下划线/连字符的字符替换为下划线（与 Data_Network 的
    prefix 过滤规则一致），阻断 `../`、`..\\`、绝对路径等路径穿越；
    再去除首尾的 `._-`，空结果退化为 "unnamed"。
    仅用于构造文件名/目录名，不要用于 API 查询参数（查询仍用原始输入）。
    """
    import re

    cleaned = re.sub(r"[^\w\-]+", "_", str(name).strip())
    return cleaned.strip("._-") or "unnamed"


def _safe_round(value, ndigits: int = 4) -> float:
    """None 安全取分：Open Targets 部分靶点的指标字段可能为 null，统一转 0.0 再四舍五入。"""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return 0.0
