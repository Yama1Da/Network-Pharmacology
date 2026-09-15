"""界面中英文文案对照表。

界面上所有文字都用 t("键名") 从这里取,如 t("app.title")。
词典本体在 translations_en.py / translations_zh.py,本文件只保留对外 API。
"""
import streamlit as st

from .translations_en import TEXT_EN
from .translations_zh import TEXT_ZH

# 默认语言:英文
DEFAULT_LANG = "en"
# 支持的语言
SUPPORTED_LANGS = {"en": "English", "zh": "中文"}

TEXT: dict[str, dict[str, str]] = {"en": TEXT_EN, "zh": TEXT_ZH}


def get_lang() -> str:
    """当前语言,默认英文。"""
    return st.session_state.get("lang", DEFAULT_LANG)


def set_lang(lang: str) -> None:
    """切换语言。"""
    if lang in TEXT:
        st.session_state["lang"] = lang


def t(key: str, **kwargs) -> str:
    """按键取当前语言文案,缺失时回退英文,再缺失返回键名。"""
    text = TEXT.get(get_lang(), {}).get(key)
    if text is None:
        text = TEXT[DEFAULT_LANG].get(key, key)
    return text.format(**kwargs) if kwargs else text
