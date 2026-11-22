# -*- coding: utf-8 -*-
"""i18n 词典导出器(薄壳):把词典合并导出为 TS 前端的 JSON。

═══ 文案修改指引 ═══
· 后端人言(错误消息/进度文案): WY/i18n/translations_{zh,en}.py —— t() 取词,双语同键
· 前端模板人言(界面导航/表格/提示/按钮): WY/i18n/translations_web_{zh,en}.py
· 改完任何词典后重跑: python -m WY.tools.export_i18n_json
· WY/web/src/locales/{zh,en}.json 是生成物(勿手改;JSON 不支持注释,指引只能写在这里)
· 进度/错误键另需在 tests/test_api.py 的 _PROG_CASES 登记渲染样例(缺键/占位符守卫)

运行方式(在仓库根目录;直接 python tools/… 会 ModuleNotFoundError):
    python -m WY.tools.export_i18n_json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from WY.i18n.en.translations_en import TEXT_EN
from WY.i18n.en.translations_web_en import TEXT_WEB_EN
from WY.i18n.zh.translations_web_zh import TEXT_WEB_ZH
from WY.i18n.zh.translations_zh import TEXT_ZH

WEB_ZH, WEB_EN = TEXT_WEB_ZH, TEXT_WEB_EN  # 旧名别名:键集断言沿用历史称呼

WEB_I18N = {"zh": WEB_ZH, "en": WEB_EN}

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "web", "src", "locales")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    # zh/en 键集合必须一致,防止某语言漏译时静默回退。
    diff = set(WEB_ZH) ^ set(WEB_EN)
    if diff:
        raise SystemExit(f"web 词典 zh/en 键不一致: {sorted(diff)}")
    for name, text in (("zh", TEXT_ZH), ("en", TEXT_EN)):
        merged = {**text, **WEB_I18N[name]}
        out = os.path.join(OUT_DIR, f"{name}.json")
        with open(out, "w", encoding="utf-8") as f:
            json.dump(merged, f, ensure_ascii=False, indent=1, sort_keys=True)
        print(f"[i18n] {out}: {len(merged)} keys")


if __name__ == "__main__":
    main()
