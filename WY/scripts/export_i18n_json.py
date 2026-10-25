# -*- coding: utf-8 -*-
"""i18n 词典导出:把 core 的两本词典(translations_zh/en)导出为 TS 前端的 JSON。

WY/web/src/locales/{zh,en}.json 由本脚本生成(勿手改,重跑即刷新);
同时并入 TS 前端自己的 web.* 命名空间(导航等界面骨架文案),
前端 vue-i18n 以这两份 JSON 为消息源。运行方式(在 WY 目录下):
    python scripts/export_i18n_json.py
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from WY.i18n.translations_en import TEXT_EN
from WY.i18n.translations_zh import TEXT_ZH

# TS 前端骨架文案(web.* 命名空间):导航与通用按钮。
# 视图内文案的键化迁移逐页进行;此处先覆盖骨架。
WEB_KEYS = {
    "web": {
        "nav": {
            "home": "主页", "analysis": "分析", "docking": "对接", "data": "数据",
            "logs": "日志", "drug_sources": "药物", "disease_sources": "疾病",
            "ppi": "PPI", "enrichment": "富集", "stars": "收藏",
            "ml": "ML", "ml_results": "ML 结果", "settings": "设置",
        },
        "common": {
            "refresh": "刷新", "all": "全选", "unselect_all": "取消全选",
            "download_selected": "下载选中 (ZIP)", "delete_selected": "删除选中",
            "selected": "已选 {n} 项", "confirm_delete": "确认删除 {n} 项?",
            "search_ph": "搜索…", "language": "Language",
        },
    }
}

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "web", "src", "locales")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, text in (("zh", TEXT_ZH), ("en", TEXT_EN)):
        web = WEB_KEYS["web"]
        merged = {**text, **{f"web.{k}.{k2}": v
                             for k, sub in web.items() for k2, v in sub.items()}}
        out = os.path.join(OUT_DIR, f"{name}.json")
        with open(out, "w", encoding="utf-8") as f:
            json.dump(merged, f, ensure_ascii=False, indent=1, sort_keys=True)
        print(f"{out}: {len(merged)} keys")


if __name__ == "__main__":
    main()
