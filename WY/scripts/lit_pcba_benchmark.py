#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""LIT-PCBA 基准验证脚本。
本脚本与 build_interactome.py为计划长期保留的两个 scripts/ 脚本。
"""

import argparse
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import config
import core.docking as dock

# ---------------- 指标计算(纯函数,可直接被 pytest 复用) ----------------


def auc_mw(labels: np.ndarray, values: np.ndarray) -> float:
    """ROC-AUC(Mann-Whitney 秩和,平均秩处理并列)。values 越大代表打分越好。"""
    labels = np.asarray(labels, dtype=int)
    values = np.asarray(values, dtype=float)
    n_a = int((labels == 1).sum())
    n_i = int((labels == 0).sum())
    if n_a == 0 or n_i == 0:
        raise ValueError("AUC 需要同时存在活性与无活性分子")
    from scipy.stats import rankdata

    ranks = rankdata(values, method="average")  # 秩 1 = 最小值(最差)
    r_a = ranks[labels == 1].sum()
    return float((r_a - n_a * (n_a + 1) / 2.0) / (n_a * n_i))


def auc_best_first(labels: np.ndarray, vina_scores: np.ndarray) -> tuple[float, float]:
    """用文档 §6 的 Excel 口径重算 AUC:按 vina 分数升序(越负越好)排名,秩 1 = 最好。
    返回 (AUC, 活性分子秩和 R)。与 auc_mw 互为独立实现,自检时互验。"""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(vina_scores, dtype=float)
    n_a = int((labels == 1).sum())
    n_i = len(labels) - n_a
    order = np.argsort(scores, kind="mergesort")  # 升序 = 最好在前
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0  # 并列取平均名次(RANK.AVG 语义)
        i = j + 1
    r_sum = float(ranks[labels == 1].sum())
    auc = (n_a * n_i + n_a * (n_a + 1) / 2.0 - r_sum) / (n_a * n_i)
    return float(auc), r_sum


def ef_at(labels: np.ndarray, vina_scores: np.ndarray, frac: float = 0.01) -> tuple[float, int]:
    """富集因子 EF@frac:前 frac 比例内的活性密度 ÷ 全库活性密度。
    k = int(frac*N)(向下取整,不足 1 取 1);并列按 (分数, CID) 确定性排序。返回 (EF, k)。"""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(vina_scores, dtype=float)
    n = len(labels)
    n_a = int((labels == 1).sum())
    if n_a == 0 or n == n_a:
        raise ValueError("EF 需要同时存在活性与无活性分子")
    k = max(1, int(frac * n))
    order = np.lexsort((np.arange(n), scores))  # 分数升序,并列按原始顺序(= CID 顺序)
    hits = int(labels[order[:k]].sum())
    return float((hits / k) / (n_a / n)), k


def bedroc(labels: np.ndarray, values: np.ndarray, alpha: float = 20.0) -> float:
    """BEDROC(alpha=20,Truchon & Bayly 2007 标准参数),用 RDKit 自带实现。
    注意:RDKit 接口为 CalcBEDROC(scores, col, alpha),col 是"真标签"所在列,
    分数本身由列表排序位置决定(必须最好在前);并列时不做同分校正,按传入顺序(= CID 序)。"""
    from rdkit.ML.Scoring import Scoring

    rows = [(float(v), int(l)) for l, v in zip(labels, values)]
    rows.sort(key=lambda r: r[0], reverse=True)  # 最好在前;RDKit 只看位置不看分数值
    return float(Scoring.CalcBEDROC(rows, 1, alpha))  # col=1 = 标签列


def permutation_auc(labels: np.ndarray, values: np.ndarray, n_perm: int = 20,
                    seed: int = 0) -> tuple[float, float]:
    """置换对照:打乱打分与标签的对应关系后 AUC 应退化到约 0.5。返回 (均值, 标准差)。"""
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_perm):
        aucs.append(auc_mw(labels, rng.permutation(values)))
    return float(np.mean(aucs)), float(np.std(aucs))


# ---------------- 数据读取 ----------------


def load_ligands(target_dir: str, split: str) -> list[dict]:
    """读取靶点目录下的活性/无活性 SMILES 清单。split=both 合并 T/V 两份(基准默认口径)。"""
    if split == "both":
        patterns = ("active*", "inactive*")
    else:
        patterns = (f"active_{split}*", f"inactive_{split}*")
    ligands: list[dict] = []
    seen: dict[str, str] = {}
    for pattern in patterns:
        matches = sorted(p for p in os.listdir(target_dir)
                         if p.startswith(pattern.rstrip("*")) and p.endswith(".smi"))
        if not matches:
            raise FileNotFoundError(f"{target_dir} 下找不到 {pattern}*.smi")
        for fname in matches:
            label = 1 if fname.startswith("active") else 0
            with open(os.path.join(target_dir, fname), encoding="utf-8") as f:
                for line in f:
                    parts = line.split()
                    if not parts or parts[0].startswith("#"):
                        continue
                    smiles, cid = parts[0], parts[1] if len(parts) > 1 else ""
                    if cid in seen:
                        if seen[cid] != smiles:
                            print(f"[警告] CID {cid} 重复且 SMILES 不一致,保留首个", file=sys.stderr)
                        continue
                    seen[cid] = smiles
                    ligands.append({"cid": cid, "smiles": smiles, "label": label,
                                    "name": f"{'act' if label else 'inact'}_{cid}"})
    if not ligands:
        raise RuntimeError(f"{target_dir} 未读到任何分子")
    return ligands


def list_templates(target_dir: str) -> list[str]:
    """列出候选受体模板(PDB 号,来自 *_protein.mol2 文件名)。"""
    return sorted(p[:-len("_protein.mol2")] for p in os.listdir(target_dir)
                  if p.endswith("_protein.mol2"))


# ---------------- 平台同款准备链 ----------------


def prepare_receptor(target_dir: str, pdb_id: str, work_dir: str) -> tuple[str, dict]:
    """mol2 -> pdb(obabel) -> pdbqt(平台 prepare_receptor_pdbqt,-xr 刚性受体)。"""
    rec_dir = os.path.join(work_dir, "receptor")
    os.makedirs(rec_dir, exist_ok=True)
    pdb_path = os.path.join(rec_dir, f"{pdb_id}.pdb")
    if not (os.path.isfile(pdb_path) and os.path.getsize(pdb_path) > 0):
        dock._run_obabel(os.path.join(target_dir, f"{pdb_id}_protein.mol2"), pdb_path)
    rec_pdbqt = dock.prepare_receptor_pdbqt(pdb_path, os.path.join(rec_dir, f"{pdb_id}_rec.pdbqt"))
    lig_pdb = os.path.join(rec_dir, f"{pdb_id}_lig.pdb")
    if not (os.path.isfile(lig_pdb) and os.path.getsize(lig_pdb) > 0):
        dock._run_obabel(os.path.join(target_dir, f"{pdb_id}_ligand.mol2"), lig_pdb)
    box = dock.infer_box_from_cofactor(lig_pdb, padding=BOX_PADDING[0])
    box_json = os.path.join(rec_dir, f"{pdb_id}_box.json")
    with open(box_json, "w", encoding="utf-8") as f:
        json.dump(box, f, ensure_ascii=False, indent=1)
    return rec_pdbqt, box


BOX_PADDING = [8.0]  # 平台 infer_box_from_cofactor 默认外扩;由 --padding 覆写(模块级便于单测打桩)


def prepare_ligand(cid: str, smiles: str, lig_dir: str) -> str:
    """SMILES -> 3D SDF(RDKit,固定种子) -> pdbqt(Meeko),与网页端完全同链;结果按 CID 缓存。"""
    os.makedirs(lig_dir, exist_ok=True)
    pdbqt_path = os.path.join(lig_dir, f"{cid}.pdbqt")
    if os.path.isfile(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
        return pdbqt_path
    sdf_path = os.path.join(lig_dir, f"{cid}.sdf")
    if not (os.path.isfile(sdf_path) and os.path.getsize(sdf_path) > 0):
        dock.smiles_to_sdf_3d(smiles, sdf_path, seed=0, n_confs=20)
    return dock.prepare_ligand_pdbqt(sdf_path, pdbqt_path)


def parse_vina_score(out_pdbqt: str) -> float | None:
    """从 vina 输出 pdbqt 的首条 REMARK VINA RESULT 取最优打分(kcal/mol,越负越好)。"""
    try:
        with open(out_pdbqt, encoding="utf-8") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    return float(line.split()[3])
    except OSError:
        pass
    return None


# ---------------- 单靶点主流程 ----------------


def run_target(target: str, args) -> dict:
    target_dir = os.path.join(args.data, target)
    if not os.path.isdir(target_dir):
        raise FileNotFoundError(f"靶点目录不存在:{target_dir}")
    templates = list_templates(target_dir)
    if not templates:
        raise FileNotFoundError(f"{target_dir} 下没有 *_protein.mol2 受体模板")
    pdb_id = args.pdb
    if pdb_id is None:
        if len(templates) > 1:
            raise SystemExit(
                f"{target} 有 {len(templates)} 个候选受体模板:{', '.join(templates)}。\n"
                "对表前请先按所对论文选定模板(指南 §5.3 第②条),用 --pdb <PDB号> 指定。")
        pdb_id = templates[0]
    if pdb_id not in templates:
        raise SystemExit(f"模板 {pdb_id} 不在 {target_dir} 中。候选:{', '.join(templates)}")

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    # 输出文件名带上模板与划分:不同受体模板/划分的结果各存各的,既不互相覆盖、也绝不允许混算
    tag = f"{target}_{pdb_id}_{args.split}"
    work_dir = os.path.join(out_dir, "work", target, pdb_id)
    lig_dir = os.path.join(out_dir, "work", target, "lig")  # 配体准备与模板无关,按靶点共享
    scores_csv = os.path.join(out_dir, f"scores_{tag}.csv")
    cfg_sig = (f"pdb={pdb_id},split={args.split},padding={BOX_PADDING[0]},"
               f"exhaustiveness={config.EXHAUSTIVENESS},num_modes={config.NUM_MODES},"
               f"seed={config.DOCKING_SEED},cpu={config.DOCKING_CPU}")

    ligands = load_ligands(target_dir, args.split)
    if args.limit:
        actives = [x for x in ligands if x["label"] == 1][:args.limit]
        inacts = [x for x in ligands if x["label"] == 0][:args.limit]
        ligands = actives + inacts

    done: dict[str, dict] = {}
    cfg_path = os.path.join(out_dir, f"scores_{tag}.cfg")
    if os.path.isfile(scores_csv):
        if not os.path.isfile(cfg_path):
            raise SystemExit(f"缺少守卫文件 {cfg_path}(scores 与 .cfg 必须成对):"
                             "请删除旧 scores 文件后重跑。")
        with open(cfg_path, encoding="utf-8") as f:
            cfg_line = f.read().strip()
        if cfg_line != f"# config: {cfg_sig}":
            raise SystemExit(
                f"{scores_csv} 由另一组参数生成（{cfg_line}），与本次"
                f"（# config: {cfg_sig}）不一致,不能混在同一文件续跑:"
                "请换 --out 输出目录,或删除旧结果文件后重跑。")
        with open(scores_csv, encoding="utf-8-sig", newline="") as f:
            lines = f.readlines()
        data = [ln for ln in lines if ln.strip() and not ln.startswith("cid,")]
        for row in csv.DictReader(data, fieldnames=["cid", "label", "score", "seconds", "status"]):
            done[row["cid"]] = row
        print(f"[续跑] {scores_csv} 已有 {len(done)} 条结果,跳过这些分子")
    # 断点续跑时重试此前非 ok 的行(vina 超时/被 Ctrl+C 杀死等瞬态失败;
    # RDKit 坏 SMILES 会再失败一次,每次开销不足一秒)
    pending = [x for x in ligands
               if x["cid"] not in done or done[x["cid"]]["status"] != "ok"]
    print(f"[{target}] 模板 {pdb_id} | 活性 {sum(1 for x in ligands if x['label'] == 1)} + "
          f"无活性 {sum(1 for x in ligands if x['label'] == 0)} | 待跑 {len(pending)} 个")

    rec_pdbqt, box = prepare_receptor(target_dir, pdb_id, work_dir)
    print(f"[{target}] 口袋 center=({box['center'][0]:.1f},{box['center'][1]:.1f},"
          f"{box['center'][2]:.1f}) size=({box['size'][0]:.1f},{box['size'][1]:.1f},"
          f"{box['size'][2]:.1f}) | 受体 {rec_pdbqt}")

    csv_exists = os.path.isfile(scores_csv)
    scores_file = open(scores_csv, "a", encoding="utf-8", newline="")
    writer = csv.DictWriter(scores_file, fieldnames=["cid", "label", "score", "seconds", "status"])
    if not csv_exists:
        # 纯标准 CSV(仅 BOM + 表头 + 数据):Excel/PyCharm/pandas 都按逗号正确解析,
        # 无指令行/注释行;参数指纹写入同名 .cfg 守卫文件,防止不同参数混入续跑
        scores_file.write("\ufeff")  # BOM:Excel 识别 UTF-8
        writer.writeheader()
        with open(cfg_path, "w", encoding="utf-8") as f:
            f.write(f"# config: {cfg_sig}\n")

    rows: list[dict] = list(done.values())
    t0 = time.time()
    try:
        for idx, lig in enumerate(pending, 1):
            row = {"cid": lig["cid"], "label": lig["label"], "score": "", "seconds": "", "status": "ok"}
            try:
                started = time.time()
                lig_pdbqt = prepare_ligand(lig["cid"], lig["smiles"], lig_dir)
                out_pdbqt = os.path.join(work_dir, "dock", f"{lig['cid']}.pdbqt")
                os.makedirs(os.path.dirname(out_pdbqt), exist_ok=True)
                dock.run_vina(rec_pdbqt, lig_pdbqt, box, out_pdbqt)
                score = parse_vina_score(out_pdbqt)
                if score is None:
                    raise RuntimeError("输出中未找到 REMARK VINA RESULT")
                row["score"], row["seconds"] = f"{score:.3f}", f"{time.time() - started:.1f}"
            except KeyboardInterrupt:
                raise
            except dock.VinaInterrupted as exc:
                # Ctrl+C 杀死了正在跑的 vina:按中断处理,当前分子不写行(续跑时自动重做)
                raise KeyboardInterrupt() from exc
            except Exception as exc:  # 单分子失败不中断整批:记录后继续
                row["status"] = f"failed: {str(exc)[:120]}"
                print(f"  [失败] CID {lig['cid']}: {str(exc)[:120]}", file=sys.stderr)
            writer.writerow(row)
            scores_file.flush()
            rows.append(row)
            ok = [r for r in rows if r["status"] == "ok" and r["score"] != ""]
            if row["status"] == "ok":
                elapsed = time.time() - t0
                eta_min = elapsed / len(ok) * (len(pending) - idx) / 60 if ok else 0
                print(f"[{target} {idx}/{len(pending)}] CID {lig['cid']} "
                      f"score={row['score']} ({row['seconds']}s) ETA {eta_min:.0f}min")
    except KeyboardInterrupt:
        scores_file.close()
        print(f"\n[中断] {target} 已完成 {len(rows)} 条,已落盘;重跑同一命令自动续跑剩余分子")
        raise  # 交由 main 立即终止整个程序,不再继续下一个靶点
    scores_file.close()

    return evaluate_and_report(target, pdb_id, tag, args, ligands, rows, box, out_dir,
                               time.time() - t0)


def evaluate_and_report(target: str, pdb_id: str, tag: str, args, ligands, rows, box,
                        out_dir, elapsed_s: float) -> dict:
    """算指标、写 metrics json 与 ROC 图、打印摘要。失败分子计入对应计数,不进指标。"""
    latest: dict[str, dict] = {}
    for r in rows:
        latest[r["cid"]] = r  # 同 cid 多行(失败重试)以最后一行为准
    ok_rows = [r for r in latest.values() if r["status"] == "ok"]
    prep_or_dock_fail = len(latest) - len(ok_rows)
    labels = np.array([int(r["label"]) for r in ok_rows], dtype=int)
    vina = np.array([float(r["score"]) for r in ok_rows], dtype=float)
    n_a = int((labels == 1).sum())
    n_i = len(labels) - n_a
    metrics: dict = {
        "target": target, "pdb_template": pdb_id, "data_dir": args.data,
        "split": args.split, "box": {"center": list(box["center"]), "size": list(box["size"]),
                                     "padding": BOX_PADDING[0]},
        "params": {"exhaustiveness": config.EXHAUSTIVENESS, "num_modes": config.NUM_MODES,
                   "cpu": config.DOCKING_CPU, "seed": config.DOCKING_SEED,
                   "timeout_s": config.DOCKING_TIMEOUT},
        "counts": {"loaded": len(ligands), "ok": len(ok_rows),
                   "failed": prep_or_dock_fail,
                   "actives_ok": n_a, "inactives_ok": n_i,
                   "note": "官方数 = loaded;failed 主要是 RDKit 无法嵌入 3D 的坏 SMILES(永久性),"
                           "vina 瞬态失败(超时/被杀)在续跑时自动重试并以最后一行为准"},
        "elapsed_s": round(elapsed_s, 1),
    }
    if n_a == 0 or n_i == 0:
        metrics["error"] = "可用分子全为同一类,无法计算指标"
    else:
        values = -vina  # 越负越好 -> 越大越好
        auc_excel, r_sum = auc_best_first(labels, vina)
        perm_mean, perm_sd = permutation_auc(labels, values, n_perm=args.permute_n)
        ef1, k1 = ef_at(labels, vina, 0.01)
        metrics["metrics"] = {
            "roc_auc": round(auc_mw(labels, values), 4),
            "roc_auc_excel_variant": round(auc_excel, 4),
            "ef1pct": round(ef1, 4), "top1pct_count": k1,
            "bedroc_alpha20": round(bedroc(labels, values, 20.0), 4),
            "permutation_auc": {"mean": round(perm_mean, 4), "sd": round(perm_sd, 4),
                                "n": args.permute_n, "seed": 0,
                                "expect": "约 0.5;明显高于则指标实现有误"},
        }
        metrics["excel_check"] = {
            "active_rank_sum_best_first": round(r_sum, 1), "n_actives": n_a, "n_inactives": n_i,
            "formula": "AUC = (n_a*n_i + n_a*(n_a+1)/2 - R) / (n_a*n_i),秩按分数升序(越负越好),并列取平均名次",
        }
        write_roc_png(os.path.join(out_dir, f"roc_{tag}.png"), labels, values, target)

    metrics_path = os.path.join(out_dir, f"metrics_{tag}.json")
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=1)
    print(f"[{target}] 指标与计数已写入 {metrics_path}")
    if "metrics" in metrics:
        m = metrics["metrics"]
        print(f"[{target}] ROC-AUC={m['roc_auc']} EF1%={m['ef1pct']} "
              f"BEDROC={m['bedroc_alpha20']} | 置换对照 AUC={m['permutation_auc']['mean']}")
    print(f"[{target}] 对表数值:官方活性 {metrics['counts']['loaded'] and sum(1 for x in ligands if x['label'] == 1)}"
          f"/成功 {n_a},无活性 {sum(1 for x in ligands if x['label'] == 0)}/成功 {n_i},失败 {prep_or_dock_fail}")
    return metrics


def write_roc_png(path: str, labels: np.ndarray, values: np.ndarray, target: str) -> None:
    """ROC 曲线图(论文图2素材);坐标轴用英文避免中文字体缺字。"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = np.argsort(-values)
    y = labels[order]
    tps = np.cumsum(y)
    fps = np.cumsum(1 - y)
    tpr = np.concatenate(([0.0], tps / max(tps[-1], 1)))
    fpr = np.concatenate(([0.0], fps / max(fps[-1], 1)))
    fig, ax = plt.subplots(figsize=(4.2, 4.0), dpi=150)
    ax.plot(fpr, tpr, color="#1f6fb4", lw=1.6, drawstyle="default")
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, ls="--")
    auc_val = auc_mw(labels, values)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"LIT-PCBA {target} (AUC={auc_val:.3f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------- 汇总与自检 ----------------


def summarize(out_dir: str) -> None:
    """扫描 out 目录下全部 metrics_*.json,打印跨靶点汇总表并写 summary.csv。"""
    paths = sorted(p for p in os.listdir(out_dir) if p.startswith("metrics_") and p.endswith(".json"))
    if not paths:
        print("没有可汇总的 metrics_*.json")
        return
    summary_path = os.path.join(out_dir, "summary.csv")
    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["target", "pdb", "split", "roc_auc", "ef1pct", "bedroc_alpha20",
                         "perm_auc_mean", "actives_ok", "inactives_ok", "failed"])
        print(f"\n{'靶点':<10} {'PDB':<6} {'划分':<5} {'AUC':>6} {'EF1%':>7} {'BEDROC':>7} "
              f"{'置换AUC':>8} {'活性':>5} {'无活性':>6} {'失败':>4}")
        for p in paths:
            with open(os.path.join(out_dir, p), encoding="utf-8") as fh:
                m = json.load(fh)
            mm = m.get("metrics", {})
            c = m.get("counts", {})
            writer.writerow([m["target"], m["pdb_template"], m.get("split", ""),
                             mm.get("roc_auc", ""),
                             mm.get("ef1pct", ""), mm.get("bedroc_alpha20", ""),
                             mm.get("permutation_auc", {}).get("mean", ""),
                             c.get("actives_ok", ""), c.get("inactives_ok", ""),
                             c.get("failed", "")])
            print(f"{m['target']:<10} {m['pdb_template']:<6} {m.get('split','-'):<5} "
                  f"{mm.get('roc_auc','-'):>6} "
                  f"{mm.get('ef1pct','-'):>7} {mm.get('bedroc_alpha20','-'):>7} "
                  f"{mm.get('permutation_auc',{}).get('mean','-'):>8} "
                  f"{c.get('actives_ok','') :>5} {c.get('inactives_ok',''):>6} {c.get('failed',''):>4}")
    print(f"汇总已写入 {summary_path}")


def selftest() -> None:
    """指标实现对已知答案的小样本断言;纯本地秒回。全部通过打印 SELFTEST PASSED。"""
    # 手算小例(文档 §6 的 5 分子表):AUC = 5/6
    labels = np.array([1, 0, 1, 0, 0])
    vina = np.array([-8.5, -7.2, -6.9, -5.5, -4.0])
    values = -vina
    assert abs(auc_mw(labels, values) - 5 / 6) < 1e-12, auc_mw(labels, values)
    auc_bf, r_sum = auc_best_first(labels, vina)
    assert abs(auc_bf - 5 / 6) < 1e-12 and abs(r_sum - (1 + 3)) < 1e-9, (auc_bf, r_sum)
    ef, k = ef_at(labels, vina, 0.01)
    assert k == 1 and abs(ef - 2.5) < 1e-12, (ef, k)  # 前 1 名是活性:(1/1)/(2/5)=2.5

    # 完美/最差分离与全并列
    n = 100
    lab = np.array([1] * 10 + [0] * 90)
    assert auc_mw(lab, np.arange(n, 0, -1.0)) == 1.0
    assert auc_mw(lab, np.arange(1, n + 1, dtype=float)) == 0.0
    assert abs(auc_mw(lab, np.full(n, 7.0)) - 0.5) < 1e-12  # 全并列 = 瞎猜
    assert auc_mw(lab, np.arange(n, 0, -1.0)) == auc_best_first(lab, -np.arange(n, 0, -1.0))[0]

    # EF1% 边界:k 向下取整不足 1 取 1;完美分离时 EF = N/n_a
    ef, k = ef_at(lab, np.arange(n, 0, -1.0) * -1.0, 0.01)  # 活性分数最负 = 最好
    assert k == 1 and abs(ef - n / 10) < 1e-9, (ef, k)
    ef, k = ef_at(lab, np.concatenate([np.full(50, -9.0), np.linspace(-8, -1, 50)]), 0.10)
    assert k == 10 and abs(ef - 10.0) < 1e-9, (ef, k)  # 前 10 名全是活性

    # BEDROC:标准化后完美分离 = 1、完全反序 = 0(RDKit 归一化口径)
    assert abs(bedroc(lab, np.arange(n, 0, -1.0)) - 1.0) < 1e-9
    assert bedroc(lab, np.arange(1, n + 1, dtype=float)) < 1e-6

    # 置换对照:打乱后应靠近 0.5
    mean, sd = permutation_auc(lab, np.arange(n, 0, -1.0), n_perm=30)
    assert 0.35 < mean < 0.65, mean

    print("SELFTEST PASSED: AUC(秩和/Excel 两种实现互验)/EF/BEDROC/置换对照 全部通过")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="LIT-PCBA 基准验证(平台同款对接链)")
    ap.add_argument("--data", help="LIT-PCBA 解压根目录(含各靶点子目录),如 D:\\lit_pcba\\AVE_unbiased")
    ap.add_argument("--target", nargs="+", help="靶点目录名,可多个(如 TP53 ESR1_ant PPARG)")
    ap.add_argument("--pdb", help="受体模板 PDB 号(多模板时必填,单模板可省)")
    ap.add_argument("--split", choices=["both", "T", "V"], default="both",
                    help="T/V 划分口径,默认 both(合并);须与对表论文一致")
    ap.add_argument("--out", default=os.path.join(
        os.environ.get("CLI_RESULTS_DIR", os.path.join("results")), "benchmark"),
        help="输出目录(默认 CLI_RESULTS_DIR/benchmark 或 ./results/benchmark)")
    ap.add_argument("--padding", type=float, default=8.0, help="口袋外扩(埃),默认与平台一致 8.0")
    ap.add_argument("--cpu", type=int, default=None,
                    help="vina 线程数;默认沿用平台设置(config.DOCKING_CPU=0,全部核心)。"
                         "笔记本过热建议 --cpu 4:温度友好,单分子约慢一倍,打分结果不受影响"
                         "(搜索深度由 exhaustiveness 决定,cpu 只影响速度);全程固定同一值即可")
    ap.add_argument("--limit", type=int, default=0, help="每类只取前 N 个分子(冒烟用)")
    ap.add_argument("--permute-n", type=int, default=20, help="置换对照次数")
    ap.add_argument("--selftest", action="store_true", help="只跑指标自检(纯本地)")
    ap.add_argument("--summarize", action="store_true", help="只汇总 --out 下已有结果")
    args = ap.parse_args(argv)

    BOX_PADDING[0] = args.padding
    if args.cpu is not None:
        config.DOCKING_CPU = args.cpu  # 只改 vina 线程数(速度/发热);打分分布由 exhaustiveness 决定
    if args.selftest:
        selftest()
        return 0
    if args.summarize:
        summarize(args.out)
        return 0
    if not args.data or not args.target:
        ap.error("--data 与 --target 必填(或用 --selftest / --summarize)")

    overall = 0
    try:
        for target in args.target:
            try:
                run_target(target, args)
            except SystemExit as exc:  # 参数/模板类守卫:跳过该靶点,继续其余靶点
                print(str(exc), file=sys.stderr)
                overall = 1
            except Exception as exc:
                print(f"[{target}] 运行失败:{exc}", file=sys.stderr)
                overall = 1
    except KeyboardInterrupt:
        print("\n[中断] 已停止全部靶点;已完成的结果均已落盘,重跑同一命令自动续跑")
        return 130
    if overall == 0:
        summarize(args.out)
    else:
        print("部分靶点失败;其余靶点的结果仍已写入输出目录", file=sys.stderr)
    return overall


if __name__ == "__main__":
    sys.exit(main())