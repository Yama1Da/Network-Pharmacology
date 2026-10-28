"""core.geo 的离线解析测试:FTP 分桶规则、series matrix 自解析、gzip 校验、下载重试兜底。

阶段三第 5 批新写。不发真实网络(下载打桩),series matrix 用合成 gz 喂解析器。
"""
import gzip
import io
import os

import pytest

from WY import config
import WY.core.geo as geo


def test_geo_bucket_rules():
    """FTP 千位分桶:去掉编号末三位拼 nnn;不足四位不补零。"""
    assert geo._geo_bucket("GSE63067") == "GSE63nnn"
    assert geo._geo_bucket("GSE193066") == "GSE193nnn"
    assert geo._geo_bucket("GPL570") == "GPLnnn"   # 恰好三位 → 全部替换
    assert geo._geo_bucket("GPL96") == "GPLnnn"    # 不足三位
    assert geo._geo_bucket("gse63067") == "GSE63nnn"  # 大小写归一


def test_series_matrix_url():
    """series matrix 直链按分桶规律拼出。"""
    assert geo.series_matrix_url("GSE63067") == (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63067/matrix/"
        "GSE63067_series_matrix.txt.gz")


def _synthetic_matrix_text(with_table: bool = True) -> str:
    """合成一份带引号 TSV 的 series matrix 文本(与 GEO 导出格式同构)。"""
    lines = [
        '!Series_title\t"Synthetic title"',
        '!Series_platform_id\t"GPL570"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"',
        '!Sample_characteristics_ch1\t"disease: case"\t"disease: case"\t"disease: control"',
        '!Sample_characteristics_ch1\t"tissue: blood"\t"tissue: blood"\t"tissue: blood"',
        '!Sample_platform_id\t"GPL570"\t"GPL570"\t"GPL570"',
    ]
    if with_table:
        lines += [
            "!series_matrix_table_begin",
            'ID_REF\t"GSM1"\t"GSM2"\t"GSM3"',
            'GENE1\t"1.0"\t"2.0"\t"3.0"',
            'GENE2\t"4.0"\t"5.0"\t"6.0"',
            "!series_matrix_table_end",
        ]
    return "\n".join(lines) + "\n"


def _write_gz(path: str, text: str) -> str:
    with open(path, "wb") as f:
        f.write(gzip.compress(text.encode("utf-8")))
    return path


def test_parse_series_matrix_full(tmp_path):
    """解析产物:样本注释表(键值对拆列)+ 数值化数据表 + 元数据 + 平台号。"""
    path = _write_gz(str(tmp_path / "GSEXXXX_series_matrix.txt.gz"),
                     _synthetic_matrix_text())
    ds = geo.parse_series_matrix(path)
    assert list(ds.samples.index) == ["GSM1", "GSM2", "GSM3"]
    # characteristics 的「键: 值」拆成独立列;单值列(tissue)也保留
    assert list(ds.samples["disease"]) == ["case", "case", "control"]
    assert set(ds.samples["tissue"]) == {"blood"}
    assert ds.table is not None and ds.table.shape == (2, 3)
    assert ds.table.loc["GENE2", "GSM3"] == 6.0
    assert ds.meta.get("Series_title") == "Synthetic title"
    assert ds.gpl == "GPL570"
    assert ds.issues == []


def test_parse_series_matrix_empty_table_recorded(tmp_path):
    """元数据型提交(无数据表):不清真地返回空表,而是记入 issues 交由调用方报错。"""
    path = _write_gz(str(tmp_path / "GSEYYYY_series_matrix.txt.gz"),
                     _synthetic_matrix_text(with_table=False))
    ds = geo.parse_series_matrix(path)
    assert ds.table.empty
    assert any("数据表为空" in i for i in ds.issues), ds.issues


def test_gzip_intact_rejects_bad_header(tmp_path):
    "不是 gzip 的文件:BadGzipFile(属 OSError)被正确接住,按不完整处理。"
    bad = tmp_path / "bad.gz"
    bad.write_bytes(b"plain text not gzip at all")
    assert geo._gzip_intact(str(bad)) is False


def test_gzip_intact_rejects_corrupt_stream(tmp_path):
    "损坏/截断的 gzip 数据流必须按不完整处理(返回 False),而不是异常炸穿调用方。"
    path = _write_gz(str(tmp_path / "ok.gz"), "content")
    raw = open(path, "rb").read()
    bad = tmp_path / "corrupt.gz"
    bad.write_bytes(raw[: len(raw) // 2])
    assert geo._gzip_intact(str(bad)) is False


def test_download_success_returns_true(tmp_path, monkeypatch):
    """下载成功(分段单次通过,内部打桩):返回 True,产物完好。"""
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRIES", 1)
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRY_DELAY", 0)

    def fake_pass(url, dest, cfg):
        _write_gz(dest, "fresh data")

    monkeypatch.setattr(geo, "_download_one_pass", fake_pass)
    dest = str(tmp_path / "x.txt.gz")
    assert geo.download_geo_file("https://x/y.gz", dest) is True
    assert geo._gzip_intact(dest)


def test_download_exhausts_raises_with_options(tmp_path, monkeypatch):
    """重试耗尽且无旧缓存:响亮报错并给出三条出路。"""
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRIES", 1)
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRY_DELAY", 0)
    monkeypatch.setattr(geo, "_download_one_pass",
                        lambda url, dest, cfg: (_ for _ in ()).throw(
                            ConnectionError("down")))
    dest = str(tmp_path / "x.txt.gz")
    with pytest.raises(RuntimeError, match="可选出路"):
        geo.download_geo_file("https://x/y.gz", dest)


def test_download_degrades_to_old_cache_returns_false(tmp_path, monkeypatch):
    """重试耗尽但本地有完好旧缓存:降级沿用,返回 False(调用方必须向用户暴露)。"""
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRIES", 1)
    monkeypatch.setattr(config, "GEO_DOWNLOAD_RETRY_DELAY", 0)
    monkeypatch.setattr(geo, "_download_one_pass",
                        lambda url, dest, cfg: (_ for _ in ()).throw(
                            ConnectionError("down")))
    dest = _write_gz(str(tmp_path / "cached.txt.gz"), "old but good")
    assert geo.download_geo_file("https://x/y.gz", dest) is False
    assert geo._gzip_intact(dest)
