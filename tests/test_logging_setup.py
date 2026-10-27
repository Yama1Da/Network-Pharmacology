"""core.logging_setup 的测试:日志初始化(控制台+按天文件)、当前日志防误删、
运行中动态调级。

阶段三第 1 批新写(原自检脚本未覆盖此模块)。LOG_DIR 由 conftest 沙箱接管。
"""
import logging
import os

from WY import config
from WY.core import logging_setup


def _snapshot_root_logger():
    """保存根 logger 的处理器与级别(测试后恢复,不污染 pytest 自身的日志设施)。"""
    root = logging.getLogger()
    return root.handlers[:], root.level


def _restore_root_logger(saved):
    handlers, level = saved
    root = logging.getLogger()
    for h in root.handlers[:]:
        root.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    for h in handlers:
        root.addHandler(h)
    root.setLevel(level)


def test_setup_logging_adds_console_and_daily_file():
    """首次初始化:挂上控制台与按天文件两个处理器,文件落在沙箱日志目录;
    重复调用不重复设置(幂等)。"""
    saved = _snapshot_root_logger()
    try:
        root = logging.getLogger()
        for h in saved[0]:
            root.removeHandler(h)
        assert not root.hasHandlers()
        logging_setup.setup_logging()
        kinds = [type(h).__name__ for h in root.handlers]
        assert "StreamHandler" in kinds and "_DailyRotatingFileHandler" in kinds, kinds
        fh = next(h for h in root.handlers
                  if type(h).__name__ == "_DailyRotatingFileHandler")
        assert os.path.dirname(fh.baseFilename) == config.LOG_DIR
        # 幂等:已有时不再追加
        logging_setup.setup_logging()
        assert len(root.handlers) == 2
    finally:
        _restore_root_logger(saved)


def test_daily_handler_never_lists_current_log(tmp_path):
    """过期日志清单(Python 3.13 语义:getFilesToDelete 返回超出 backupCount 的
    待删最旧文件):两个轮转旧文件 + backupCount=1 → 只列更旧的;清单里绝不含
    当前正在写入的 app.log(误删当前日志在 Windows 上会直接报错)。"""
    current = tmp_path / "app.log"
    current.write_text("today", encoding="utf-8")
    (tmp_path / "app.log.2026-01-01").write_text("newer", encoding="utf-8")
    (tmp_path / "app.log.2025-01-01").write_text("older", encoding="utf-8")
    handler = logging_setup._DailyRotatingFileHandler(
        str(current), when="midnight", backupCount=1, encoding="utf-8")
    try:
        leftovers = [os.path.basename(p) for p in handler.getFilesToDelete()]
    finally:
        handler.close()
    assert leftovers == ["app.log.2025-01-01"], leftovers
    assert "app.log" not in leftovers


def test_set_logging_level_runtime_change():
    """动态调级:根 logger 与全部处理器同步生效;非法级别名原样不动。"""
    root = logging.getLogger()
    saved_level = root.level
    saved_handler_levels = [h.level for h in root.handlers]
    try:
        logging_setup.set_logging_level("DEBUG")
        assert root.level == logging.DEBUG
        assert all(h.level == logging.DEBUG for h in root.handlers)
        logging_setup.set_logging_level("not-a-level")  # 非法名:不改任何东西
        assert root.level == logging.DEBUG
    finally:
        root.setLevel(saved_level)
        for h, lv in zip(root.handlers, saved_handler_levels):
            h.setLevel(lv)
