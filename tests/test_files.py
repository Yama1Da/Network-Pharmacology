"""core.files 的测试:缓存/下载/日志三类列表的过滤规则,四种 zip 打包。

阶段三第 1 批新写(原自检脚本未覆盖此模块)。目录路径全部走 conftest 沙箱。
"""
import os
import zipfile

from WY import config
from WY.core.files import (
    list_cache_files,
    list_download_dirs,
    list_log_files,
    make_dir_zip,
    make_dirs_zip,
    make_fig_zip,
    make_ppi_zip,
)


def _mkfile(path: str, content: str = "x") -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return path


def test_list_cache_files_filters(tmp_path):
    """缓存列表:只列文件不列目录、跳过隐藏文件与 .bak/.tmp、可按扩展名过滤。"""
    category = os.path.join(config.CACHE_DIR, "testcat")
    _mkfile(os.path.join(category, "a.csv"))
    _mkfile(os.path.join(category, "b.csv.bak"))
    _mkfile(os.path.join(category, ".hidden.csv"))
    os.makedirs(os.path.join(category, "subdir"))
    out = list_cache_files("testcat")
    assert list(out["File name"]) == ["a.csv"], list(out["File name"])
    # 扩展名过滤
    _mkfile(os.path.join(category, "c.json"))
    out_json = list_cache_files("testcat", exts=(".json",))
    assert list(out_json["File name"]) == ["c.json"]
    # 目录不存在 → 空表
    assert list_cache_files("no_such_cat").empty


def test_list_download_dirs(tmp_path):
    """下载列表:只列子目录不列散文件。"""
    os.makedirs(os.path.join(config.DOWNLOADS_DIR, "aspirin_cancer_data"))
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "loose.txt"))
    out = list_download_dirs()
    assert list(out["Directory"]) == ["aspirin_cancer_data"]


def test_list_log_files(tmp_path):
    """日志列表:app.log 与 app.log.日期 形态都列出,其余文件不列。"""
    _mkfile(os.path.join(config.LOG_DIR, "app.log"))
    _mkfile(os.path.join(config.LOG_DIR, "app.log.2026-01-01"))
    _mkfile(os.path.join(config.LOG_DIR, "other.txt"))
    out = list_log_files()
    assert set(out["File name"]) == {"app.log", "app.log.2026-01-01"}


def test_make_ppi_zip(tmp_path):
    """结果打包:存在的文件按 basename 入包;不存在与 directory 键跳过。"""
    f1 = _mkfile(str(tmp_path / "ppi_edges.csv"), "a,b\n1,2\n")
    f2 = _mkfile(str(tmp_path / "summary.json"), "{}")
    data = make_ppi_zip({"directory": str(tmp_path), "edges": f1,
                         "summary": f2, "missing": str(tmp_path / "nope.csv")})
    import io

    zf = zipfile.ZipFile(io.BytesIO(data))
    assert set(zf.namelist()) == {"ppi_edges.csv", "summary.json"}


def test_make_fig_zip(tmp_path):
    """一图一钮打包:PNG+SVG 双格式;缺一个就只装另一个。"""
    png = _mkfile(str(tmp_path / "fig.png"), "png")
    svg = _mkfile(str(tmp_path / "fig.svg"), "svg")
    zf = zipfile.ZipFile(io_bytes(make_fig_zip(png, svg)))
    assert set(zf.namelist()) == {"fig.png", "fig.svg"}
    zf2 = zipfile.ZipFile(io_bytes(make_fig_zip(png, str(tmp_path / "nope.svg"))))
    assert zf2.namelist() == ["fig.png"]


def io_bytes(data: bytes):
    import io

    return io.BytesIO(data)


def test_make_dir_zip(tmp_path):
    """整目录打包:目录内全部文件入包;目录不存在 → 空 zip。"""
    _mkfile(str(tmp_path / "d1" / "a.txt"))
    _mkfile(str(tmp_path / "d1" / "b.txt"))
    zf = zipfile.ZipFile(io_bytes(make_dir_zip(str(tmp_path / "d1"))))
    assert set(zf.namelist()) == {"a.txt", "b.txt"}
    zf_empty = zipfile.ZipFile(io_bytes(make_dir_zip(str(tmp_path / "no_dir"))))
    assert zf_empty.namelist() == []


def test_make_dirs_zip(tmp_path):
    """批量目录打包:每个目录作为 zip 内一层文件夹,保留目录名前缀;缺失目录跳过。"""
    _mkfile(str(tmp_path / "run1" / "a.csv"))
    _mkfile(str(tmp_path / "run2" / "sub" / "b.csv"))
    zf = zipfile.ZipFile(io_bytes(make_dirs_zip(
        {"run1": str(tmp_path / "run1"), "run2": str(tmp_path / "run2"),
         "nope": str(tmp_path / "no_dir")})))
    # zipfile 入包时把 os.sep 归一为 "/",弧名断言用正斜杠
    names = set(zf.namelist())
    assert names == {"run1/a.csv", "run2/sub/b.csv"}
