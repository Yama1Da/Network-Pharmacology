"""本机工具链探测：AutoDock Vina 与 Open Babel 的可执行文件/目录发现逻辑。

从 config.py 拆出的零依赖叶子（只 import 标准库）：config 在导入时调用本模块取值
填 VINA_DIR/VINA_EXE 等常量；config.apply_api_env 在用户经齿轮面板改目录后调用
find_exe 重探对应可执行文件。探测结果只反映"当前环境"：换机器后改 API.env 里的
VINA_DIR/OBABEL_DIR 即可，无需改代码。"""

import os


def find_exe(tool_dir: str, prefix: str, exclude: str = "") -> str:
    """目录里按文件名前缀找可执行文件(不限定平台与版本文件名)。

    多个匹配时取字典序最后一个;exclude 用于主程序排除 vina_split;找不到返回空串。
    (第一参数在旧实现里叫 vina_dir,拆出后更名 tool_dir:它泛指"该工具所在目录",
    对 Vina 与 OBABEL 通用。)"""
    if not tool_dir or not os.path.isdir(tool_dir):
        return ""
    names = sorted(
        f for f in os.listdir(tool_dir)
        if f.lower().startswith(prefix)
        and not (exclude and f.lower().startswith(exclude))
        and os.path.isfile(os.path.join(tool_dir, f))
    )
    return os.path.join(tool_dir, names[-1]) if names else ""


def find_vina_dir() -> str:
    """探测 Vina 程序所在目录(跨平台):环境变量 VINA_DIR 优先,其次 PATH。

    VINA_DIR 或 PATH 里填目录即可(也可填可执行文件全路径)。"""
    env_dir = os.environ.get("VINA_DIR")
    if env_dir:
        if os.path.isfile(env_dir):                                   # 填的是可执行文件全路径
            return os.path.dirname(env_dir)
        if find_exe(env_dir, "vina", exclude="vina_split"):           # 填的是目录
            return env_dir

    # vina 所在目录已加入 PATH(按 vina* 前缀匹配,兼容各平台发布文件名)
    for path_dir in os.environ.get("PATH", "").split(os.pathsep):
        if find_exe(path_dir.strip().strip('"'), "vina", exclude="vina_split"):
            return path_dir.strip().strip('"')

    return ""  # 找不到时,运行对接会报错提示设置 VINA_DIR


def find_obabel_dir() -> str:
    """探测 Open Babel(obabel)所在目录(跨平台):环境变量 OBABEL_DIR 优先,其次 PATH。

    OBABEL_DIR 或 PATH 里填目录即可(也可填可执行文件全路径)。"""
    env_dir = os.environ.get("OBABEL_DIR")
    if env_dir:
        if os.path.isfile(env_dir):                                   # 填的是可执行文件全路径
            return os.path.dirname(env_dir)
        if find_exe(env_dir, "obabel"):                               # 填的是目录
            return env_dir

    # obabel 所在目录已加入 PATH
    for path_dir in os.environ.get("PATH", "").split(os.pathsep):
        if find_exe(path_dir.strip().strip('"'), "obabel"):
            return path_dir.strip().strip('"')

    return ""  # 找不到时,运行转换会报错提示设置 OBABEL_DIR
