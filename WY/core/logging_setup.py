"""日志设置:运行记录同时显示在屏幕上,并按天存入日志文件(保留两年)。"""

import logging
import logging.handlers
import os


class _DailyRotatingFileHandler(logging.handlers.TimedRotatingFileHandler):
    """按天存日志的文件(每天零点切换到新文件)。
    修正过期日志的自动清理:只删除真正过期的旧日志文件,不误删当前正在写入的日志(否则在 Windows 上会报错)。
    """

    def getFilesToDelete(self) -> list[str]:
        files = super().getFilesToDelete()
        current = os.path.basename(self.baseFilename)
        return [f for f in files if os.path.basename(f) != current]


# DEBUG 模式的常驻噪音源,永久钳到 WARNING:matplotlib 开 DEBUG 会逐个枚举系统字体
# (一开几百行),PIL 吐图像解码细节,joblib 吐并行批处理细节,均与排查本项目问题无关。
# 子级 logger(如 matplotlib.font_manager)默认 NOTSET,钳父级即一并覆盖。
# 网络层有意放行(不在此列):urllib3/requests_cache 的连接与缓存日志本身就在 DEBUG 级,
# 不会污染 INFO 视图,而是 DEBUG 按钮的核心价值——还原下载/远程抓取的现场。
_NOISY_LOGGERS = ("matplotlib", "PIL", "joblib")


def _clamp_noisy_loggers() -> None:
    """把第三方噪音 logger 钳到 WARNING;只设子 logger 级别,不挂处理器,随日志级别开关保持不变。"""
    for name in _NOISY_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


def setup_logging() -> None:
    """初始化日志:同时输出到屏幕和按天滚动的日志文件,保留两年;重复调用不会重复设置。
    过期日志自动删除;只记录重要信息,屏蔽画图库的琐碎调试信息(_NOISY_LOGGERS);日志目录创建失败时只在屏幕输出,不影响程序运行。
    """
    from WY import config

    logger = logging.getLogger()
    if logger.hasHandlers():
        return

    level = getattr(logging, str(config.LOG_LEVEL).upper(), logging.INFO)
    if not isinstance(level, int):
        level = logging.INFO
    logger.setLevel(level)
    _clamp_noisy_loggers()
    formatter = logging.Formatter(config.LOG_FORMAT, config.LOG_DATE_FORMAT)

    # 把 warnings.warn 的输出并入日志框架:否则第三方库警告(sklearn/scipy 等)
    # 只裸打到 stderr,无时间戳/级别前缀,app.log 也没有,事后无法排查。
    # 桥接后走 py.warnings 记录器,WARNING 级,格式与本项目日志一致。
    logging.captureWarnings(True)

    # 屏幕输出(方便实时查看运行情况)
    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # 文件输出(按天一个文件,保留两年后自动删除)
    try:
        os.makedirs(config.LOG_DIR, exist_ok=True)
        log_file = os.path.join(config.LOG_DIR, "app.log")
        file_handler = _DailyRotatingFileHandler(
            log_file,
            when="midnight",                              # 每天零点切文件
            backupCount=config.LOG_RETENTION_DAYS,        # 保留两年，超期自动删除
            encoding="utf-8",
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except OSError:
        # 日志文件创建失败:改为只在屏幕输出,不影响程序运行
        logger.warning("日志文件创建失败，降级为仅控制台输出: %s", config.LOG_DIR)


def set_logging_level(level_name: str) -> None:
    """切换运行中的根日志级别,使后续日志立即生效。"""
    level = getattr(logging, str(level_name).upper(), None)
    if not isinstance(level, int):
        return
    logger = logging.getLogger()
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)
    logging.info("Log level changed to %s", logging.getLevelName(level))
