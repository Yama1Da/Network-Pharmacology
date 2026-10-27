"""core.management 的测试:UI 设置落盘白名单、API.env 设置链路;env_discovery 工具链探测。

迁移自 selfcheck_docking.py(18 项中的 3 项)。ui_settings.json 路径由 conftest
沙箱接管;os.environ 操纵改用 monkeypatch.setenv/delenv(测试结束自动还原)。
"""
import json
import os

from WY import config
from WY.core.management import save_ui_settings


def test_ui_settings_roundtrip():
    """口袋半径作为全局分析参数:在齿轮白名单内,经 save_ui_settings 落盘到
    ui_settings.json,且多次写入按键合并(对接页与齿轮面板互不覆盖)。"""
    assert "POCKET_CUTOFF" in config.UI_SETTINGS_KEYS, config.UI_SETTINGS_KEYS
    save_ui_settings({"POCKET_CUTOFF": 4.5})
    save_ui_settings({"DGIDB_MIN_SCORE": 0.1})
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        saved = json.load(f)
    assert saved == {"POCKET_CUTOFF": 4.5, "DGIDB_MIN_SCORE": 0.1}, saved
    # 白名单外的键(如黑名单以外的任意字符串值)不落盘
    save_ui_settings({"NOT_IN_WHITELIST": 1.0})
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        assert "NOT_IN_WHITELIST" not in json.load(f)


def test_api_env_roundtrip(monkeypatch, tmp_path):
    """API.env 设置链路(save_api_env):白名单过滤、值消毒(换行/井号拒收)、set_key 落盘
    保留注释与其他行、运行时回填含派生值联动(NCBI key 生效 → E-UTILS 节流间隔切 0.10)。"""
    from WY.core.management import save_api_env

    base = str(tmp_path)
    env_path = os.path.join(base, "API.env")
    with open(env_path, "w", encoding="utf-8") as f:
        f.write("# 注释行必须保留\nCLI_RESULTS_DIR=D:/old\n")
    orig = {k: getattr(config, k) for k in config.API_ENV_KEYS}
    orig_interval = config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"]
    try:
        monkeypatch.setattr(config, "BASE_DIR", base)
        applied = save_api_env({
            "DisGeNET_API_KEY": "  abc123  ",   # 消毒:去首尾空白后生效
            "NCBI_API_KEY": "bad\nvalue",       # 消毒:含换行拒收
            "VINA_DIR": "",                     # 空值跳过(留空=不修改)
            "NOT_IN_WHITELIST": "x",            # 白名单外丢弃
            "OBABEL_DIR": '"D:/obabel dir"',    # 消毒:剥掉成对首尾引号(复制路径常带引号)
            "CLI_RESULTS_DIR": "D:/new",        # 已有行被覆盖
        })
        assert applied == ["DisGeNET_API_KEY", "OBABEL_DIR", "CLI_RESULTS_DIR"], applied
        assert config.DisGeNET_API_KEY == "abc123" and config.CLI_RESULTS_DIR == "D:/new"
        assert config.OBABEL_DIR == "D:/obabel dir", config.OBABEL_DIR
        with open(env_path, encoding="utf-8") as f:
            text = f.read()
        assert "# 注释行必须保留" in text
        assert "abc123" in text and "D:/new" in text and "D:/old" not in text
        assert "NOT_IN_WHITELIST" not in text and "bad" not in text
        # 引号已剥离:文件里不应出现带引号的值
        assert '"D:/obabel dir"' not in text
        # 运行时回填联动派生值:NCBI key 生效 → E-UTILS 节流间隔切 0.10
        save_api_env({"NCBI_API_KEY": "ncbi_key_456"})
        assert config.NCBI_API_KEY == "ncbi_key_456"
        assert config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] == 0.10
    finally:
        # save_api_env 会原地修改运行时 config(API_ENV_KEYS 各值与节流间隔),
        # monkeypatch 管不到,按自检脚本同款手法显式恢复
        config.apply_api_env(orig)
        config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] = orig_interval


def test_env_discovery(monkeypatch, tmp_path):
    """工具链探测(env_discovery,自 config 拆出):find_exe 的前缀匹配/exclude/非文件与
    空目录安全空串,find_vina_dir 的 环境变量>PATH 回退与"可执行文件全路径"识别、
    都找不到时返回空串。打桩 PATH 只含测试目录,不受本机真实安装干扰。"""
    from WY import env_discovery

    base = str(tmp_path)
    vina_dir = os.path.join(base, "vina_bin")
    os.makedirs(vina_dir)
    for name in ("vina_1.2.7_win.exe", "vina_split_1.2.7_win.exe",
                 "obabel.exe", "readme.txt"):
        open(os.path.join(vina_dir, name), "wb").close()
    # find_exe:前缀匹配 + exclude 排除 vina_split + 多匹配取字典序最后一个
    assert env_discovery.find_exe(vina_dir, "vina", exclude="vina_split") == \
        os.path.join(vina_dir, "vina_1.2.7_win.exe")
    assert env_discovery.find_exe(vina_dir, "vina_split") == \
        os.path.join(vina_dir, "vina_split_1.2.7_win.exe")
    assert env_discovery.find_exe(vina_dir, "obabel") == \
        os.path.join(vina_dir, "obabel.exe")
    # 目录本身(其名以 vina 开头但是目录)、空参、不存在的目录:安全返回空串
    assert env_discovery.find_exe(base, "vina") == ""
    assert env_discovery.find_exe("", "vina") == ""
    assert env_discovery.find_exe(os.path.join(base, "nope"), "vina") == ""
    # find_vina_dir:环境变量优先;填可执行文件全路径时返回其所在目录
    monkeypatch.setenv("VINA_DIR", vina_dir)
    assert env_discovery.find_vina_dir() == vina_dir
    monkeypatch.setenv("VINA_DIR", os.path.join(vina_dir, "vina_1.2.7_win.exe"))
    assert env_discovery.find_vina_dir() == vina_dir
    # PATH 回退:打桩期间 PATH 只含该目录,防本机真实安装干扰
    monkeypatch.delenv("VINA_DIR", raising=False)
    monkeypatch.setenv("PATH", vina_dir)
    assert env_discovery.find_vina_dir() == vina_dir
    # 都找不到:返回空串(运行对接时报错提示设置 VINA_DIR)
    monkeypatch.setenv("PATH", base)
    assert env_discovery.find_vina_dir() == ""
    # find_obabel_dir 同规则:环境变量指向目录即命中
    monkeypatch.setenv("OBABEL_DIR", vina_dir)
    assert env_discovery.find_obabel_dir() == vina_dir
