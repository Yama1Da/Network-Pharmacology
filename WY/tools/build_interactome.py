#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""命令行脚本:从 STRING 数据库下载并构建全蛋白背景互作网络的本地缓存。

用法（在 WY 目录下）：
    python tools/build_interactome.py                 # 构建 Human（默认参数）
    python tools/build_interactome.py --force         # 已构建也强制重建
    python tools/build_interactome.py --confidence 700
    python tools/build_interactome.py --keep-raw      # 保留原始下载文件（调试用）

产物为背景网络的边表(基因对 + 互作得分)、基因连接数表和构建记录。
同参数重跑会覆盖原文件;构建成功后默认删除原始下载文件。
"""

import argparse
import os
import sys
import time


from WY.core.interactome import (  # noqa: E402
    build_interactome,
    is_built,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build STRING interactome local cache (M2)")
    parser.add_argument("--species", default="Human",
                        help="物种（默认 Human；Mouse/Rat 需先在 config.INTERACTOME_ENABLED_SPECIES 放开）")
    parser.add_argument("--confidence", type=int, default=None,
                        help="STRING combined_score 阈值（0-1000），默认取 config.INTERACTOME_CONFIDENCE")
    parser.add_argument("--force", action="store_true",
                        help="已构建时也强制重建（默认跳过并提示）")
    parser.add_argument("--keep-raw", action="store_true",
                        help="保留下载的原始 gz 文件（调试用，默认构建成功后删除）")
    args = parser.parse_args()

    if is_built(args.species, args.confidence) and not args.force:
        print(f"[skip] interactome for {args.species} already built "
              f"(use --force to rebuild).")
        return 0

    t0 = time.time()
    print(f"[build] species={args.species} confidence={args.confidence or 'config default'} ...")
    meta = build_interactome(
        species=args.species,
        confidence=args.confidence,
        keep_raw=args.keep_raw,
    )
    dt = time.time() - t0
    from WY import config

    print("[done] nodes=%d edges=%d version=%s confidence>=%d elapsed=%.1fs"
          % (meta["n_nodes"], meta["n_edges"], meta["string_version"],
             meta["confidence"], dt))
    print("[done] cache dir:", config.INTERACTOME_CACHE_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
