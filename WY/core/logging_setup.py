"""日志系统初始化：控制台 + 按天轮转文件（保留两年），幂等。"""

import logging
import logging.handlers
import os


class _DailyRotatingFileHandler(logging.handlers.TimedRotatingFileHandler):
    """按天轮转的文件 handler（每天零点切文件）。

    继承 TimedRotatingFileHandler 并修复其清理逻辑：原版 getFilesToDelete()
    会把"当前正在写入的文件"（如 app.log）也计入待删除列表，在 Windows 上
    删除被占用的文件会抛 PermissionError，导致两年后的过期清理失效。
    这里在删除前排除当前文件，使 backupCount 清理只删真正的旧轮转文件。
    """

    def getFilesToDelete(self) -> list[str]:
        files = super().getFilesToDelete()
        current = os.path.basename(self.baseFilename)
        return [f for f in files if os.path.basename(f) != current]


def setup_logging() -> None:
    """初始化日志系统：控制台 + 每日文件，保留两年，幂等。

    - 日志目录 config.LOG_DIR（WY/Static/log）
    - 【按天轮转】TimedRotatingFileHandler 每天零点自动把日志写入新的
      轮转文件 app.log.YYYY-MM-DD（服务器不重启也能正确切分到当天文件）
    - 【保留两年】backupCount = LOG_RETENTION_DAYS（730），超过自动删除
      最旧的轮转文件，取代手写的 _cleanup_old_logs
    - 【轻量版】root logger 级别设为 INFO：屏蔽 matplotlib/PIL 等第三方库的
      DEBUG 噪音（字体查找、图片流），仅录项目自身 INFO 及以上日志；
      若后续需要精细控制（项目 DEBUG + 第三方 WARNING），再升级为
      "项目专用 logger（logging.getLogger('WY')）+ root 降级" 的方案
    - 幂等：Streamlit 每次交互都会重跑脚本，用 hasHandlers 守卫避免重复挂 handler
    - 目录创建失败时降级为仅控制台输出，不抛异常
    """
    import config

    logger = logging.getLogger()
    if logger.hasHandlers():
        return  # 已初始化过

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(config.LOG_FORMAT, config.LOG_DATE_FORMAT)

    # 控制台 handler（实时排障）
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # 文件 handler（按天轮转 + 保留两年自动清理）
    try:
        os.makedirs(config.LOG_DIR, exist_ok=True)
        log_file = os.path.join(config.LOG_DIR, "app.log")
        file_handler = _DailyRotatingFileHandler(
            log_file,
            when="midnight",                              # 每天零点切文件
            backupCount=config.LOG_RETENTION_DAYS,        # 保留两年，超期自动删除
            encoding="utf-8",
        )
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except OSError:
        # 目录/文件创建失败：降级为仅控制台，不影响应用
        logger.warning("日志文件创建失败，降级为仅控制台输出: %s", config.LOG_DIR)
