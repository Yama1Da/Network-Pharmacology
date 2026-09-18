"""通用文本 / 数值小工具。"""

import os
import threading

from .caching import load_json_file, save_json_file

# safe 名 → 原始名 映射文件的进程内互斥锁(并发抓取时防止同时读写弄坏文件)
_name_map_lock = threading.Lock()


# Windows 保留设备名:查询名恰为这些词时直接建文件必失败,前缀加下划线规避。
# 映射表(remember_name_mapping)仍按 safe 名记录原始名,管理页还原不受影响。
_WINDOWS_RESERVED = frozenset(
    ["CON", "PRN", "AUX", "NUL"]
    + [f"COM{i}" for i in range(1, 10)]
    + [f"LPT{i}" for i in range(1, 10)]
)


def safe_name(name: str) -> str:
    """把用户输入(药物名/疾病名/基因等)转成能安全用作文件名的形式:文件名不允许的符号统一替换成下划线,
    防止输入内容被当作路径、破坏文件组织;转换结果为空时用 "unnamed";Windows 保留设备名加下划线前缀。
    只用于文件命名,查询时仍使用原始输入。
    """
    import re

    cleaned = re.sub(r"[^\w\-]+", "_", str(name).strip()).strip("._-") or "unnamed"
    if cleaned.upper() in _WINDOWS_RESERVED:
        cleaned = f"_{cleaned}"
    return cleaned


def _name_map_file() -> str:
    """safe 名 → 原始名 映射文件路径。"""
    import config

    return os.path.join(config.CACHE_DIR, "name_map.json")


def remember_name_mapping(original: str) -> None:
    """记录 safe_name(名) → 原始名 的映射,供管理页按文件名还原原始查询名。
    safe_name 是有损转换,管理页按文件名反推名字重新抓取时查不到数据,故各抓取函数存文件时顺带记录映射;无损转换不写盘。
    """
    import config

    original = str(original).strip()
    key = safe_name(original)
    if key == original:
        return
    with _name_map_lock:  # 读-改-写整段互斥,避免并发时丢映射
        path = _name_map_file()
        mapping = load_json_file(path, {}) or {}
        if mapping.get(key) == original:
            return
        mapping[key] = original
        save_json_file(path, mapping, indent=1)


def load_name_map() -> dict[str, str]:
    """读取完整的 safe 名 → 原始名 映射;文件缺失或损坏时返回空映射。"""
    return load_json_file(_name_map_file(), {}) or {}


def resolve_original_name(stem: str) -> str | None:
    """按 safe 名还原记录的原始名;没有映射记录(旧文件/无损名字)时返回 None。"""
    original = load_name_map().get(stem)
    return str(original) if original else None


def safe_round(value, ndigits: int = 4) -> float:
    """安全取数:部分靶点在数据库中没有指标数值,统一按 0 处理后再四舍五入。"""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return 0.0
