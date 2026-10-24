#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""界面小组件自检脚本:列表分页状态机(page_slice)与翻页条(render_pager)。

纯本地运行,不发网络请求。page_slice 只读写 st.session_state、render_pager 只
创建控件,故用最小 FakeSt 打桩隔离 bare-mode 差异,断言状态与边界逻辑本身。
运行方式(在 WY 目录下):
    python scripts/selfcheck_ui.py
"""

import os
import sys


from WY.scripts.selfcheck_common import patch_attr, run

import WY.ui.widgets as widgets


class _Col:
    """列容器:与 streamlit 一致,既是 with 上下文,自身也直接可调 caption/button。"""

    def __init__(self, st: "_FakeSt"):
        self._st = st

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def caption(self, text, *args, **kwargs):
        self._st.captions.append(str(text))

    def button(self, *args, **kwargs):
        self._st.buttons.append(kwargs)
        return False


class _FakeSt:
    """page_slice/render_pager 用到的 streamlit 面最小化:会话字典 + 控件记录。

    button 记录 kwargs(断言 disabled/key)且默认返回 False(未点击);
    rerun 记录调用次数(真实实现抛 RerunException,这里只验证"会触发重绘")。
    """

    def __init__(self):
        self.session_state = {}
        self.buttons: list[dict] = []
        self.captions: list[str] = []
        self.rerun_calls = 0

    def columns(self, spec):
        return [_Col(self) for _ in spec]

    def button(self, *args, **kwargs):
        # 行内控件在 with 列容器块里仍经 st.button 调用(容器只决定渲染位置)
        self.buttons.append(kwargs)
        return False

    def rerun(self):
        self.rerun_calls += 1


def t_page_first_page():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        result = widgets.page_slice("k", 25, query="a")
    assert result == (0, 10, 0, 3), result  # 25 条 → 3 页,首页取前 10


def t_page_single_page_counts():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        assert widgets.page_slice("k", 10, query="")[3] == 1  # 恰好一页
        assert widgets.page_slice("k2", 1, query="")[3] == 1


def t_page_empty_total_still_one_page():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        result = widgets.page_slice("k", 0, query="")
    assert result[2:] == (0, 1), result  # 空列表给出合法页码,空态由调用方渲染


def t_page_query_change_resets():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        widgets.page_slice("k", 25, query="a")
        st.session_state["k_page"] = 2
        assert widgets.page_slice("k", 25, query="a")[2] == 2  # 同查询保持页码
        assert widgets.page_slice("k", 25, query="b")[2] == 0  # 新查询从头浏览


def t_page_clamp_when_shrinks():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        widgets.page_slice("k", 25, query="a")
        st.session_state["k_page"] = 2
        result = widgets.page_slice("k", 5, query="a")  # 删除/筛选后末页越界
    assert result == (0, 10, 0, 1), result


def t_pager_hidden_single_page():
    st = _FakeSt()  # 未实现 columns:单页若误渲染翻页条会立刻 AttributeError
    with patch_attr(widgets, "st", st):
        widgets.render_pager("k", page=0, n_pages=1)


def t_pager_layout_and_bounds():
    st = _FakeSt()
    with patch_attr(widgets, "st", st):
        widgets.render_pager("k", page=0, n_pages=3)
    assert len(st.buttons) == 2, st.buttons
    prev, nxt = st.buttons
    assert (prev["key"], nxt["key"]) == ("k_page_prev", "k_page_next")
    assert prev["disabled"] is True and nxt["disabled"] is False  # 首页时上一页禁用
    assert any("1 / 3" in c for c in st.captions), st.captions


def t_pager_next_click_mutates_page():
    st = _FakeSt()

    def click(*args, **kwargs):
        return kwargs.get("key") == "k_page_next"  # 模拟点击"下一页"

    st.button = click
    with patch_attr(widgets, "st", st):
        widgets.render_pager("k", page=0, n_pages=3)
    assert st.session_state.get("k_page") == 1 and st.rerun_calls == 1


def t_i18n_page_keys_complete():
    from WY.i18n import TEXT, t

    for key in ("common.page_prev", "common.page_next", "common.page_status"):
        for lang, table in TEXT.items():
            assert key in table, (lang, key)
    # 文案占位符可用(默认语言为 en,缺失键会回退成键名本身)
    assert t("common.page_status", cur=1, total=3) == "Page 1 / 3"


def main() -> None:
    tests = [
        ("page_slice first page & page count", t_page_first_page),
        ("page_slice single-page counts", t_page_single_page_counts),
        ("page_slice empty total still one page", t_page_empty_total_still_one_page),
        ("page_slice query change resets page", t_page_query_change_resets),
        ("page_slice clamp when list shrinks", t_page_clamp_when_shrinks),
        ("render_pager hidden on single page", t_pager_hidden_single_page),
        ("render_pager layout & bounds", t_pager_layout_and_bounds),
        ("render_pager next click mutates page", t_pager_next_click_mutates_page),
        ("i18n page keys complete", t_i18n_page_keys_complete),
    ]
    run(tests, "ui widgets selfcheck")


if __name__ == "__main__":
    main()
