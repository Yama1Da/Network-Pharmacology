#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""架构自检脚本:分层约定与导入卫生的静态检查(AST 扫描,不发网络请求)。

守护的约定(见 core/__init__.py docstring 与本文件内各检查 docstring):
    段 1(导入统一前):
      1. core 零 streamlit —— core/、config.py、env_discovery.py、i18n/ 不得 import streamlit
      2. 下层不引上层 —— 上述下层文件不得引用 pages/ui/App/cli/scripts
      3. sys.path 补丁白名单 —— sys.path.insert/append 只允许出现在白名单入口文件
    段 2(已激活):
      4. sys.path 零容忍(白名单收紧为本脚本自身的条件 bootstrap)
      5. 单身份断言 —— 导入 WY.cli 全链路后,sys.modules 不得出现裸名身份的
         config/core/ui/i18n 模块(同一文件以两个名字加载两遍会让模块级锁/缓存翻倍)

运行方式(在 WY 目录下):
    python scripts/selfcheck_architecture.py
"""

import ast
import os
import sys

# 条件 bootstrap:可编辑安装(npap)已提供 WY 包可见性时零 sys.path 操作;
# 未安装环境直接以脚本方式运行本文件时,才把项目根放进 sys.path。
import importlib.util as _ilu
if _ilu.find_spec("WY") is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 项目根

from WY.scripts.selfcheck_common import run

# ======== 段 2 已激活:导入统一完成(sys.path 补丁清零 + 全仓 WY. 绝对导入) ========
STRICT = True

# sys.path 白名单(段 2):仅剩两处"条件 bootstrap"——未安装环境下的兜底,
# 正常可编辑安装环境里它们的条件都不命中,等效于零 sys.path 操作。
SYS_PATH_WHITELIST = {
    "WY/App.py",
    "WY/scripts/selfcheck_architecture.py",
}

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 项目根
WY_DIR = os.path.join(ROOT, "WY")

# 下层文件范围:业务/数据/基础层,全部不得触碰 streamlit 与上层
LOWER_PREFIXES = ("WY/core/", "WY/i18n/")
LOWER_FILES = {"WY/config.py", "WY/env_discovery.py", "WY/__init__.py"}

# 上层(界面/入口)的规范模块名:下层引用任何一个都算越界
UPPER = ("WY.pages", "WY.ui", "WY.App", "WY.cli", "WY.scripts")

# 裸名 → WY.* 规范名(就地运行约定下的导入身份;段 2 统一后裸名不再出现,映射保留无害)
BARE_TO_WY = {
    "config": "WY.config", "core": "WY.core", "ui": "WY.ui", "i18n": "WY.i18n",
    "pages": "WY.pages", "App": "WY.App", "cli": "WY.cli", "scripts": "WY.scripts",
    "env_discovery": "WY.env_discovery",
    "selfcheck_common": "WY.scripts.selfcheck_common",
}


def _py_files():
    """遍历 WY/ 下全部 .py(跳过 __pycache__),产出 (repo相对posix路径, 绝对路径)。"""
    for dirpath, dirnames, filenames in os.walk(WY_DIR):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            if fn.endswith(".py"):
                abs_path = os.path.join(dirpath, fn)
                rel = os.path.relpath(abs_path, ROOT).replace(os.sep, "/")
                yield rel, abs_path


def _is_lower(rel: str) -> bool:
    return rel.startswith(LOWER_PREFIXES) or rel in LOWER_FILES


def _canonical_imports(rel: str, abs_path: str):
    """解析一个文件的全部 import,产出 (顶层名, 规范全名) 列表。

    裸名按 BARE_TO_WY 映射到 WY.*;相对导入按文件所属包解析为绝对位置
    (含函数内延迟导入——它们同样是依赖)。解析失败的文件按无导入处理并回传告警。"""
    out, errors = [], []
    try:
        tree = ast.parse(open(abs_path, encoding="utf-8").read())
    except SyntaxError as e:
        return out, [f"{rel}: SyntaxError {e}"]
    parts = rel[:-3].split("/")
    pkg = parts[:-1]  # "WY/core/ppi.py"→["WY","core"];包 __init__ 同式
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                pairs, errs = _canon(a.name, 0, pkg)
                out.extend(pairs)
                errors.extend(errs)
        elif isinstance(node, ast.ImportFrom):
            pairs, errs = _canon(node.module or "", node.level, pkg)
            out.extend(pairs)
            errors.extend(errs)
    return out, errors


def _canon(target: str, level: int, pkg: list):
    """把一个 import 目标规约为 [(顶层名, 规范全名)];相对导入按包解析为 WY.* 绝对位置,
    逃出 WY 包时报错。裸名按 BARE_TO_WY 映射(段 2 统一后裸名消失,映射保留无害)。"""
    if level:
        base = pkg[:len(pkg) - (level - 1)]
        full = ".".join(base + ([target] if target else []))
        if not full.startswith("WY"):
            return [], [f"relative import escapes WY package: {full!r}"]
        return [(full.split(".")[0], full)], []
    top = target.split(".")[0]
    full = BARE_TO_WY[top] + target[len(top):] if top in BARE_TO_WY else target
    return [(full.split(".")[0], full)], []


def _lower_imports(files):
    """下层文件的全部规范导入;返回 (导入列表, 错误列表)。"""
    imports, errors = [], []
    for rel, abs_path in files:
        if _is_lower(rel):
            imps, errs = _canonical_imports(rel, abs_path)
            imports.append((rel, imps))
            errors.extend(errs)
    return imports, errors


def t_core_no_streamlit():
    """检查 1:core/、config.py、env_discovery.py、i18n/ 不得 import streamlit。
    业务/数据层不被 UI 框架污染,是 Streamlit 可整体替换为 TS+Vue 前端的前提;
    本检查守护它不被任何一次'顺手'破坏。"""
    imports, errors = _lower_imports(_py_files())
    assert not errors, "; ".join(errors)
    bad = [f"{rel}: {top}" for rel, imps in imports
           for top, full in imps
           if top == "streamlit" or full == "streamlit" or full.startswith("streamlit.")]
    assert not bad, "下层文件 import streamlit(违反分层铁律): " + "; ".join(bad)


def t_lower_not_upper():
    """检查 2:下层文件不得引用 pages/ui/App/cli/scripts。依赖方向是单行道
    (上层→下层),反向引用会把业务层反向绑定到某个界面上,是分层瓦解的第一步。"""
    imports, errors = _lower_imports(_py_files())
    assert not errors, "; ".join(errors)
    bad = [f"{rel}: {full}" for rel, imps in imports
           for top, full in imps
           for u in UPPER
           if full == u or full.startswith(u + ".")]
    assert not bad, "下层引用上层(依赖方向越界): " + "; ".join(bad)


def _sys_path_hits(abs_path: str) -> list:
    """找出文件里所有 sys.path.insert/append 调用的行号。"""
    hits = []
    tree = ast.parse(open(abs_path, encoding="utf-8").read())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in ("insert", "append")
                and isinstance(node.func.value, ast.Attribute)
                and node.func.value.attr == "path"
                and isinstance(node.func.value.value, ast.Name)
                and node.func.value.value.id == "sys"):
            hits.append(node.lineno)
    return hits


def t_sys_path_whitelist():
    """检查 3:sys.path 补丁只允许出现在白名单入口文件。白名单之外的新增意味着
    又一个'碰巧能用'的启动方式——双身份单例事故的种子。段 2 收紧为一项
    (本脚本自身的条件 bootstrap),由 STRICT=True 启用的零容忍检查接管。"""
    bad = []
    for rel, abs_path in _py_files():
        if rel in SYS_PATH_WHITELIST:
            continue
        hits = _sys_path_hits(abs_path)
        if hits:
            bad.append(f"{rel}: lines {hits}")
    assert not bad, "sys.path 补丁出现在白名单之外(改用 WY. 绝对导入): " + "; ".join(bad)


def t_sys_path_zero():
    """(段 2)检查 4:sys.path 补丁零容忍,全仓只允许本脚本自身的条件 bootstrap。"""
    allowed = {"WY/scripts/selfcheck_architecture.py"}
    bad = []
    for rel, abs_path in _py_files():
        if rel in allowed:
            continue
        hits = _sys_path_hits(abs_path)
        if hits:
            bad.append(f"{rel}: lines {hits}")
    assert not bad, "sys.path 补丁应已清零: " + "; ".join(bad)


def t_single_identity():
    """(段 2)检查 5:单身份断言。导入 WY.cli 全链路后,sys.modules 不得出现任何
    裸名身份的 config/core/ui/i18n 模块——同一文件以两个名字加载两遍,模块级锁/
    缓存/注册表静默翻倍,是本项目历史上最大的结构性隐患。"""
    import importlib
    importlib.import_module("WY.cli")
    bare = [m for m in sys.modules
            if m in ("config", "core", "ui", "i18n")
            or m.startswith(("core.", "ui.", "i18n."))]
    assert not bare, f"裸名身份模块仍在被加载: {bare}"


def t_guard_self_test():
    """检查器自测:用合成的违规/合法文件喂同一条检查管线,确认护栏自身能抓到
    streamlit 渗入、下层引上层、白名单外 sys.path,且不冤枉合法写法——
    防'护栏先失效再放行违规'。合成文件落进真实 WY/core/,让包解析走真实路径。"""
    synth_bad = ("import streamlit\n"
                 "from ui.widgets import x\n"
                 "def f():\n"
                 "    import sys\n"
                 "    sys.path.insert(0, '.')\n")
    synth_ok = ("import logging\n"
                "from i18n import t\n"
                "def f():\n"
                "    import config\n"
                "    from .caching import load\n")
    real_bad = "WY/core/_fake_bad.py"
    real_ok = "WY/core/_fake_ok.py"
    try:
        for rel, src in ((real_bad, synth_bad), (real_ok, synth_ok)):
            with open(os.path.join(ROOT, rel), "w", encoding="utf-8") as f:
                f.write(src)
        # 违规文件:三条检查各抓到各自的靶子
        imps_bad, errs_bad = _lower_imports([(real_bad, os.path.join(ROOT, real_bad))])
        assert not errs_bad, errs_bad
        pairs_bad = imps_bad[0][1]
        assert ("streamlit", "streamlit") in pairs_bad, pairs_bad
        assert any(full == "WY.ui.widgets" for _, full in pairs_bad), pairs_bad
        assert _sys_path_hits(os.path.join(ROOT, real_bad)) == [5]
        # 合法文件:core 内相对导入 + 下层引用(config/i18n/标准库)不触发任何检查
        imps_ok, errs_ok = _lower_imports([(real_ok, os.path.join(ROOT, real_ok))])
        assert not errs_ok, errs_ok
        pairs_ok = imps_ok[0][1]
        assert ("logging", "logging") in pairs_ok, pairs_ok
        assert ("WY", "WY.i18n") in pairs_ok, pairs_ok
        assert ("WY", "WY.config") in pairs_ok, pairs_ok
        assert ("WY", "WY.core.caching") in pairs_ok, pairs_ok
        assert not any(full == "streamlit" or full.startswith("streamlit.")
                       for _, full in pairs_ok), pairs_ok
        assert not any(full == u or full.startswith(u + ".")
                       for _, full in pairs_ok for u in UPPER), pairs_ok
        assert _sys_path_hits(os.path.join(ROOT, real_ok)) == []
    finally:
        for rel in (real_bad, real_ok):
            path = os.path.join(ROOT, rel)
            if os.path.exists(path):
                os.remove(path)


def main():
    tests = [
        ("core 零 streamlit", t_core_no_streamlit),
        ("下层不引上层", t_lower_not_upper),
        ("检查器自测", t_guard_self_test),
    ]
    if STRICT:
        tests += [
            ("sys.path 条件 bootstrap 白名单", t_sys_path_zero),
            ("单身份断言(导入 WY.cli)", t_single_identity),
        ]
    else:
        tests += [("sys.path 补丁白名单", t_sys_path_whitelist)]
    run(tests, "architecture selfcheck")


if __name__ == "__main__":
    main()
