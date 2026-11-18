"""后端文案对照表(core/api/CLI 生成用户可见字符串时取词)。

后端用 t("键名") 取词,如 t("dock.err_invalid_smiles");TS 前端的文案不在这里,
在 tools/export_i18n_json.py 的 web.* 命名空间(经该脚本导出为 web/src/locales/*.json)。
词典本体在 translations_en.py / translations_zh.py,本文件只保留对外 API。

本包不依赖任何 Web 框架:当前语言由宿主经 set_lang_provider 注入
读取函数,未注入时(离线脚本/线程等无界面环境)恒为默认语言,
因此 core 等模块可以安全 import 本包。
"""
from collections.abc import Callable

from .translations_en import TEXT_EN
from .translations_zh import TEXT_ZH

# 默认语言:英文
DEFAULT_LANG = "en"

TEXT: dict[str, dict[str, str]] = {"en": TEXT_EN, "zh": TEXT_ZH}

# 语言读取函数:宿主注入(如读 st.session_state 的语言开关),缺省恒为默认语言
_lang_provider: Callable[[], str] = lambda: DEFAULT_LANG


def set_lang_provider(getter: Callable[[], str]) -> None:
    """注册"当前语言"的读取函数,如 lambda: st.session_state.get("lang", DEFAULT_LANG)。"""
    global _lang_provider
    _lang_provider = getter


def get_lang() -> str:
    """当前语言,默认英文。"""
    return _lang_provider()


def t(key: str, **kwargs) -> str:
    """按键取当前语言文案,缺失时回退英文,再缺失返回键名并留痕日志
    (页面上出现裸键名时,日志里能直接定位是哪个键、哪种语言缺失)。"""
    text = TEXT.get(get_lang(), {}).get(key)
    if text is None:
        text = TEXT[DEFAULT_LANG].get(key, key)
        if text == key:
            import logging
            logging.getLogger(__name__).warning(
                "i18n 缺失键: %s (lang=%s) —— 若为新增词条,重启服务以加载新翻译表",
                key, get_lang())
    return text.format(**kwargs) if kwargs else text
