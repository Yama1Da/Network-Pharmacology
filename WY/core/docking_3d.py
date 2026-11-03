"""分子对接的交互式 3D 视图导出。

由受体 PDB 与对接构象生成自包含的 3Dmol.js 交互视图 HTML(页面预览与结果归档共用),
3Dmol.js 库文件的下载缓存与 BSD 许可文本也在此管理。
"""

import logging
import os

import requests

from .caching import save_text_atomic
from .docking_interactions import _WATER_RESNAMES, _run_obabel

# 拉取 3Dmol.js 库文件时随请求标识来源(与 docking 里 RCSB 请求的 UA 常量各自独立)
_3DMOL_FETCH_HEADERS = {"User-Agent": "NetworkPharmacologyPlatform/1.0"}


# ---------- 交互式 3D 视图导出 ----------

# 3Dmol.js 库的下载地址:官方源在前,jsdelivr 镜像兜底
_3DMOL_JS_URLS = (
    "https://3dmol.org/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol@2/build/3Dmol-min.js",
)

# 3Dmol.js 的许可文本(BSD-3-Clause,取自官方仓库 LICENSE):导出 HTML 内嵌或引用其库时
# 随文件一并提供,满足许可证"再分发须保留版权声明与许可文本"的要求。
_3DMOL_LICENSE_TEXT = """3Dmol.js incorporates code from GLmol, Three.js, and jQuery and is licensed under a BSD-3-Clause license.

Copyright (c) 2014, University of Pittsburgh
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""

_3D_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  body { margin: 0; font-family: "Segoe UI", Arial, sans-serif; }
  #cap { padding: 8px 14px; background: #f5f7fa; border-bottom: 1px solid #ddd;
         font-size: 13px; color: #333; }
  #lic { display: block; margin-top: 2px; font-size: 11px; color: #777; }
  #bar { padding: 4px 14px; border-bottom: 1px solid #eee; }
  #focus { cursor: pointer; }
  #hint { margin-left: 10px; font-size: 12px; color: #777; }
  #viewer { width: 100vw; height: calc(100vh - 88px); position: relative; }
</style>
__LIB__
</head>
<body>
<div id="cap">__CAPTION__<span id="lic">__LICENSE_NOTE__</span></div>
<div id="bar"><button id="focus">__FOCUS_BTN__</button><span id="hint">再次点击复位 · 双击视图任意处也可复位 · double-click the view to reset</span></div>
<div id="viewer"></div>
<script>
(function () {
  var receptor = __RECEPTOR__;
  var ligand = __LIGAND__;
  var het = __HET__;
  var pocket = __POCKET__;
  var interactions = __INTERACTIONS__;
  var viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
  viewer.addModel(receptor, "pdb");
  viewer.addModel(ligand, "pdb");
  viewer.setStyle({ model: 1 }, { stick: { colorscheme: "greenCarbon" } });
  if (het) {
    viewer.addModel(het, "pdb");
    viewer.setStyle({ model: 2 }, { stick: { colorscheme: "greyCarbon", opacity: 0.5, radius: 0.15 } });
  }
  // 口袋残基按链分组选择(无链号时不限定 chain)
  var byChain = {};
  (pocket || []).forEach(function (r) {
    (byChain[r.chain] = byChain[r.chain] || []).push(r.resi);
  });
  // 受体卡通样式统一由此函数设置:faded=true 时卡通整体半透明(虚化周围蛋白),
  // 突出口袋棍状、相互作用虚线与距离标注;虚化只是显示开关,信息不删减,双击复位即可还原。
  // 口袋残基的卡通与棍状在一条语句里显式声明完整样式,不依赖多次 setStyle 的合并行为
  function applyProteinStyles(faded) {
    var cartoon = faded ? { color: "spectrum", opacity: 0.3 } : { color: "spectrum" };
    viewer.setStyle({ model: 0 }, { cartoon: cartoon });
    Object.keys(byChain).forEach(function (ch) {
      var sel = { model: 0, resi: byChain[ch] };
      if (ch) { sel.chain = ch; }
      viewer.setStyle(sel, { cartoon: cartoon,
                             stick: { colorscheme: "orangeCarbon", radius: 0.2 } });
    });
    viewer.render();
  }
  applyProteinStyles(false);
  if (pocket && pocket.length) {
    pocket.forEach(function (r) {
      viewer.addLabel(r.resn + "-" + r.resi, {
        position: { x: r.pos[0], y: r.pos[1], z: r.pos[2] },
        fontSize: 11, fontColor: "black",
        showBackground: true, backgroundColor: "white",
        backgroundOpacity: 0.7
      });
    });
  }
  // 相互作用虚线圆柱:蓝=氢键 紫=盐桥 灰=疏水(用虚线圆柱而非 1px 线,白底下清晰);
  // 氢键与盐桥在中点标注距离
  var lineColors = { hbond: "#0f62fe", salt: "#c026d3", hydro: "#808080" };
  (interactions || []).forEach(function (it) {
    viewer.addCylinder({
      start: { x: it.pos1[0], y: it.pos1[1], z: it.pos1[2] },
      end: { x: it.pos2[0], y: it.pos2[1], z: it.pos2[2] },
      radius: 0.06, dashed: true, dashLength: 0.3, gapLength: 0.22,
      color: lineColors[it.kind] || "#777777"
    });
    if (it.kind !== "hydro") {
      viewer.addLabel(it.dist.toFixed(1) + "\\u00C5", {
        position: { x: (it.pos1[0] + it.pos2[0]) / 2,
                    y: (it.pos1[1] + it.pos2[1]) / 2,
                    z: (it.pos1[2] + it.pos2[2]) / 2 },
        fontSize: 10, fontColor: "#333333",
        showBackground: true, backgroundColor: "white",
        backgroundOpacity: 0.6
      });
    }
  });
  viewer.zoomTo();
  viewer.render();
  // 聚焦口袋:以配体为中心、拉远到能框住口袋残基,并把周围蛋白卡通虚化;
  // 按钮是开关:再次点击(或双击视图任意处)复位为整体视图并取消虚化
  var focusBtn = document.getElementById("focus");
  var focusLabel = focusBtn.textContent;
  var resetLabel = "复位视图 Reset view";
  var focused = false;
  function setFocus(on) {
    focused = on;
    applyProteinStyles(on);
    if (on) {
      viewer.zoomTo({ model: 1 });
      viewer.zoom(0.55, 400);
    } else {
      viewer.zoomTo();
    }
    focusBtn.textContent = on ? resetLabel : focusLabel;
    viewer.render();
  }
  focusBtn.onclick = function () { setFocus(!focused); };
  document.getElementById("viewer").ondblclick = function () { setFocus(false); };
})();
</script>
<!--
__3DMOL_LICENSE__
-->
</body>
</html>
"""


def _load_3dmol_js() -> str | None:
    """读取 3Dmol.js 库源码:本地有缓存直接用,没有则下载一次并长期缓存。
    两个下载地址都失败时返回 None,HTML 将改用在线引用(打开时需联网)。
    """
    from WY import config

    cache_path = os.path.join(config.DOCKING_CACHE_DIR, "3Dmol-min.js")
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 100_000:
        try:
            with open(cache_path, encoding="utf-8") as f:
                return f.read()
        except OSError:
            logging.exception("cached 3Dmol.js unreadable, re-downloading")
    for url in _3DMOL_JS_URLS:
        try:
            resp = requests.get(url, timeout=config.API_TIMEOUT,
                                headers=_3DMOL_FETCH_HEADERS)
            resp.raise_for_status()
            if "$3Dmol" not in resp.text or "</script" in resp.text:
                raise RuntimeError(f"unexpected content from {url}")
            os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
            # 原子写:下载中途被打断不会留下半截 JS 被上面"大小达标即信任"的判定永久复用
            save_text_atomic(cache_path, resp.text)
            logging.info("Cached 3Dmol.js (%d KB) from %s", len(resp.text) // 1024, url)
            return resp.text
        except Exception:
            logging.exception("3Dmol.js download failed: %s", url)
    return None


def pose_pdbqt_to_pdb(pose_pdbqt: str, out_path: str) -> str:
    """把对接构象从 pdbqt 转成通用 PDB 格式(供各查看器显示);失败时抛 RuntimeError。"""
    _run_obabel(pose_pdbqt, out_path)
    return out_path


def build_interactive_3d_html(receptor_pdb: str, ligand_pdb: str, caption: str = "",
                              lib_js: str | None = None,
                              pocket: list[dict] | None = None,
                              het_pdb: str | None = None,
                              interactions: list[dict] | None = None,
                              focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """生成自包含的交互式 3D 视图 HTML,样式与页面内预览一致(受体卡通彩虹色、配体绿色棍状、结合位点残基橙色棍状并标注名称、相互作用虚线、共晶配体灰色半透明参照),返回 (html 文本, 是否离线可用)。
    pocket 传入 residues_near_ligand 的结果时高亮结合位点残基;interactions 传入 interactions_near_ligand 的结果时画相互作用虚线并标注距离;
    het_pdb 传入共晶配体的 HETATM 文本时以半透明灰棍叠加作天然口袋参照;focus_label 为「聚焦口袋」按钮文案;
    lib_js 传入 3Dmol.js 源码时内嵌进文件双击即可离线打开,传 None 改用在线引用(打开时需联网)。
    """
    import html as _html
    import json

    if lib_js and "</script" not in lib_js:
        lib_tag = f"<script>\n{lib_js}\n</script>"
        offline_ok = True
    else:
        lib_tag = f'<script src="{_3DMOL_JS_URLS[0]}"></script>'
        offline_ok = False
    caption = _html.escape(caption)
    # 3Dmol.js 许可标注:内嵌时"内嵌",CDN 引用时"加载";完整许可文本以注释随文件提供
    note = ("内嵌 3Dmol.js (BSD-3-Clause)" if offline_ok
            else "本页面通过 CDN 加载 3Dmol.js (BSD-3-Clause)")
    out = (_3D_HTML_TEMPLATE
           .replace("__3DMOL_LICENSE__", _3DMOL_LICENSE_TEXT)
           .replace("__LICENSE_NOTE__",
                    f"{note} · Copyright (c) 2014, University of Pittsburgh · https://3dmol.org")
           .replace("__LIB__", lib_tag)
           .replace("__TITLE__", caption or "3D view")
           .replace("__CAPTION__", caption)
           .replace("__FOCUS_BTN__", _html.escape(focus_label))
           .replace("__INTERACTIONS__", json.dumps(interactions or []))
           .replace("__POCKET__", json.dumps(pocket or []))
           .replace("__HET__", json.dumps(het_pdb))
           .replace("__RECEPTOR__", json.dumps(receptor_pdb))
           .replace("__LIGAND__", json.dumps(ligand_pdb)))
    return out, offline_ok


def cofactor_hetatm_pdb(pdb_path: str, receptor_text: str | None = None) -> str:
    """提取受体 PDB 中的非水 HETATM(共晶配体/辅因子)行,拼成 PDB 文本,供 3D 渲染叠加天然口袋参照;
    无共晶配体时返回空字符串。receptor_text 传入已读取的 PDB 全文时不再重复读盘。"""
    if receptor_text is None:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            receptor_text = f.read()
    return "\n".join(ln for ln in receptor_text.splitlines()
                     if ln.startswith("HETATM")
                     and ln[17:20].strip() not in _WATER_RESNAMES)


def interactive_3d_html(pdb_path: str, pose_pdbqt: str, caption: str = "",
                        pocket: list[dict] | None = None,
                        receptor_text: str | None = None,
                        pose_pdb: str | None = None,
                        interactions: list[dict] | None = None,
                        focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """由受体 PDB 与所选构象 pdbqt 生成可保存的交互式 3D 视图 HTML。
    pocket 传入口袋残基列表时一并嵌入;interactions 传入相互作用列表时一并嵌入;
    受体含共晶配体(非水 HETATM)时以半透明灰棍叠加,过滤规则与页面预览一致。
    receptor_text 传入已读取的受体 PDB 全文、pose_pdb 传入已转换好的构象 PDB 路径时直接复用,
    供页面在同一次刷新里让预览与保存 HTML 共用一次转换、一次读盘;不传时保持原有的自足行为。
    """
    if pose_pdb is None:
        pose_pdb = os.path.join(os.path.dirname(pose_pdbqt), "pose_view.pdb")
        pose_pdbqt_to_pdb(pose_pdbqt, pose_pdb)
    if receptor_text is None:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            receptor_text = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=receptor_text)
    return build_interactive_3d_html(receptor_text, ligand_pdb, caption,
                                     lib_js=_load_3dmol_js(),
                                     pocket=pocket, het_pdb=het_pdb or None,
                                     interactions=interactions, focus_label=focus_label)


def save_interactive_3d_html(pdb_path: str, pose_pdbqt: str, out_path: str,
                             caption: str = "",
                             pocket: list[dict] | None = None,
                             interactions: list[dict] | None = None,
                             focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """把交互式 3D 视图保存为单个 HTML 文件(同名覆盖,便于同一构象反复保存);返回 (文件路径, 是否离线可用)。"""
    html_text, offline_ok = interactive_3d_html(pdb_path, pose_pdbqt, caption,
                                                pocket=pocket, interactions=interactions,
                                                focus_label=focus_label)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html_text)
    logging.info("Saved interactive 3D view: %s (offline=%s)", out_path, offline_ok)
    return out_path, offline_ok
