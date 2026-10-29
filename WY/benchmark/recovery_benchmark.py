# -*- coding: utf-8 -*-
"""平台级 benchmark:恢复型验证(Recovery Validation)。

设计(评估报告 §8.2):
- 金标准:Open Targets 的已知药物-疾病适应症对(APPROVAL 级);
- 防泄漏:构建疾病基因池时,剔除该药物的直接靶点证据,再跑全流程;
- 考题:平台能否把已知药物(或其靶点)重新排回 RWR/邻近度的前列;
- 指标:Top-K 回收率 + 经验 p 值(与随机排名比较)。

用法(在 WY 目录下):
    python benchmark/recovery_benchmark.py --drug metformin --disease "type 2 diabetes mellitus"
    python benchmark/recovery_benchmark.py --drug aspirin --disease "cardiovascular disease"
"""

import argparse
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import requests

import config
from i18n import t

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

_OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
_HEADERS = {"User-Agent": "Mozilla/5.0"}


def _config():
    return config


def _ot_drug_indications(chembl_id: str) -> dict:
    """查 Open Targets 获取药物的适应症(MONDO ID + 疾病名)和已知靶点。"""
    query = """
    query($chemblId: String!) {
      drug(chemblId: $chemblId) {
        id
        name
        maximumClinicalStage
        indications {
          count
          rows {
            maxClinicalStage
            disease {
              id
              name
            }
          }
        }
        linkedTargets {
          rows {
            target {
              approvedSymbol
            }
          }
        }
      }
    }
    """
    r = requests.post(_OT_URL,
                      json={"query": query, "variables": {"chemblId": chembl_id}},
                      timeout=30, headers=_HEADERS)
    r.raise_for_status()
    d = r.json().get("data", {}).get("drug", {})
    if not d:
        raise RuntimeError(f"Open Targets 未找到药物: {chembl_id}")
    inds = [{"mondo": row["disease"]["id"], "name": row["disease"]["name"],
             "stage": row.get("maxClinicalStage", "")}
            for row in d.get("indications", {}).get("rows", [])]
    targets = [row["target"]["approvedSymbol"]
               for row in d.get("linkedTargets", {}).get("rows", [])]
    return {"chembl_id": chembl_id, "name": d.get("name", ""),
            "max_phase": d.get("maximumClinicalStage", ""),
            "indications": inds, "known_targets": targets}


def _ot_disease_targets(mondo_id: str, exclude_targets: list) -> list:
    """查 Open Targets 获取疾病的关联靶点,排除指定药物的直接靶点(防泄漏)。"""
    import config
    from core.http_client import request_with_retry

    query = """
    query($eid: String!, $index: Int!, $size: Int!) {
      disease(ekgsId: $eid) {
        id
      }
      associatedTargets(diseaseId: $eid, index: $index, size: $size) {
        count
        rows {
          target {
            approvedSymbol
          }
          score
        }
      }
    }
    """
    # 使用 disease 查询
    query2 = """
    query($id: String!, $index: Int!, $size: Int!) {
      associatedTargets(diseaseId: $id, index: $index, size: $size) {
        count
        rows {
          target {
            approvedSymbol
          }
          score
        }
      }
    }
    """
    all_rows = []
    index = 0
    page_size = 100
    while True:
        r = requests.post(_OT_URL,
                          json={"query": query2,
                                "variables": {"id": mondo_id, "index": index,
                                              "size": page_size}},
                          timeout=30, headers=_HEADERS)
        r.raise_for_status()
        data = r.json().get("data", {}).get("associatedTargets", {})
        if not data.get("rows"):
            break
        for row in data["rows"]:
            symbol = row["target"]["approvedSymbol"]
            score = row.get("score", 0)
            # 防泄漏:跳过该药物的直接靶点
            if symbol not in exclude_targets:
                all_rows.append({"gene": symbol, "score": score})
        index += 1
        if index * page_size >= data.get("count", 0):
            break
    return all_rows


def recovery_ranking(ranked_genes: list, known_drug_targets: list) -> dict:
    """计算已知药物靶点在平台排名中的回收情况。"""
    rank_map = {g: i + 1 for i, g in enumerate(ranked_genes)}
    rankings = []
    for t in known_drug_targets:
        if t in rank_map:
            rankings.append({"target": t, "rank": rank_map[t], "found": True})
        else:
            rankings.append({"target": t, "rank": None, "found": False})
    found = [r for r in rankings if r["found"]]
    if not found:
        return {"found_count": 0, "top10": 0, "top20": 0, "top50": 0, "detail": rankings}
    ranks = sorted(r["rank"] for r in found)
    return {
        "found_count": len(found),
        "best_rank": ranks[0],
        "median_rank": sorted(ranks)[len(ranks)//2] if ranks else None,
        "top10": sum(1 for r in ranks if r <= 10),
        "top20": sum(1 for r in ranks if r <= 20),
        "top50": sum(1 for r in ranks if r <= 50),
        "detail": rankings,
    }


def permutation_p_value(observed_best_rank: int, n_all: int, n_known: int,
                        n_perm: int = 1000, seed: int = 0) -> float:
    """经验 p 值:随机抽 n_known 个位置的排名,统计出现 ≤ observed_best_rank 的比例。"""
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        sample = rng.choice(n_all, size=min(n_known, n_all), replace=False)
        if min(sample) + 1 <= observed_best_rank:
            count += 1
    return (count + 1) / (n_perm + 1)  # +1 校正


def main():
    ap = argparse.ArgumentParser(description="恢复型 benchmark:平台能否把已知药物重新排回来")
    ap.add_argument("--chembl", required=True, help="药物的 ChEMBL ID")
    ap.add_argument("--disease-mondo", required=True, help="疾病的 MONDO ID(如 MONDO_0005148)")
    ap.add_argument("--disease-name", required=True, help="疾病名称(用于展示)")
    ap.add_argument("--n-perm", type=int, default=1000, help="置换次数(默认 1000)")
    ap.add_argument("--out", help="结果输出目录")
    args = ap.parse_args()

    print(f"恢复型 benchmark | {args.disease_name}")
    print(f"药物 ChEMBL ID: {args.chembl}")
    print(f"疾病 MONDO ID: {args.disease_mondo}")
    print()

    # 1) 查 Open Targets:药物的适应症 + 已知靶点
    print("[1/4] 查询 Open Targets:药物适应症与已知靶点 …")
    drug_info = _ot_drug_indications(args.chembl)
    known_targets = drug_info["known_targets"]
    print(f"    药物: {drug_info['name']} | maxPhase: {drug_info['max_phase']}")
    print(f"    已知靶点: {len(known_targets)} 个")
    # 确认该药确实有这个疾病的适应症
    matched_ind = [i for i in drug_info["indications"] if args.disease_mondo in i.get("mondo", "")]
    if matched_ind:
        print(f"    ✓ 确认适应症: {matched_ind[0]['name']}")
    else:
        print(f"    ⚠ 未在 Open Targets 中确认 {args.disease_name} 适应症")

    # 2) 查疾病关联靶点(排除该药的直接靶点 = 防泄漏)
    print("[2/4] 查询疾病关联靶点(排除该药直接靶点)…")
    disease_targets = _ot_disease_targets(args.disease_mondo, known_targets)
    print(f"    排除后疾病靶点: {len(disease_targets)} 个")

    # 3) 跑平台 RWR/邻近度 → 得到排名
    # 这里用平台已有的 metformin/T2D 结果(如果药物匹配)或手动跑平台分析
    # M1 阶段为简化验证,直接用 Open Targets 的疾病靶点做 RWR 种子,
    # 并检查该药靶点在最终排名中的位置
    print("[3/4] 运行平台分析并评估回收 …")
    # 从平台 CLI 调用(或直接导入)—— 此处用 open targets 数据集做 RWR 种子
    # 实际实现:用平台已有的 RWR 排名输出(results 下)或调用 core 的 RWR 函数
    # 简化:用疾病靶点做 RWR 种子,检查药物靶点在非种子基因中的排名
    from core.targets import get_common_targets

    drug_target_set = set(known_targets)
    disease_genes = [r["gene"] for r in disease_targets]
    all_targets = list(drug_target_set) + disease_genes
    print(f"    总靶点池(药物靶点 ∪ 疾病靶点): {len(all_targets)} 个")

    # RWR 排名(简化:如果平台有输出则读取,否则用网络邻近度评估)
    # 实际实现需要调平台 RWR
    print("    (此处需要平台 RWR 输出或实时调用;详见下一版)")

    # 4) 输出
    out_dir = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "results", "recovery_benchmark")
    os.makedirs(out_dir, exist_ok=True)
    summary = {
        "module": "recovery-benchmark",
        "drug": drug_info["name"], "chembl_id": args.chembl,
        "disease": args.disease_name, "mondo_id": args.disease_mondo,
        "known_targets": known_targets, "n_known_targets": len(known_targets),
        "disease_targets_count": len(disease_targets),
        "confirm_indication": bool(matched_ind),
        "note": "M1 阶段仅完成数据验证;RWR 调用和排名计算待 M2+ 实现",
    }
    with open(os.path.join(out_dir, "recovery_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    print(f"\n结果: {out_dir}/recovery_summary.json")


if __name__ == "__main__":
    main()
