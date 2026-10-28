"""docking 抓取链与 geo 检索链补测(全部打桩,零网络)。

深水区冲线批 II。
"""
import pandas as pd
import pytest

from WY import config
import WY.core.docking as dock
import WY.core.expression as expr
import WY.core.geo as geo
import WY.core.http_client as http_client


class _FakeResp:
    def __init__(self, text="", status=200):
        self.status_code = status
        self.text = text

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"http {self.status_code}")


# ---------- docking:配体链 ----------


def test_resolve_cid_digit_passthrough():
    """输入本身是编号:直接返回整数,零网络。"""
    assert dock.resolve_cid("4091") == 4091


def test_resolve_cid_not_found_vs_unavailable(monkeypatch):
    """查无化合物(404)→ None;服务端 503 重试耗尽 → PubChemUnavailable。"""
    import pubchempy
    monkeypatch.setattr(config, "PUBCHEM_BUSY_RETRIES", 0)
    monkeypatch.setattr(config, "PUBCHEM_BUSY_BACKOFF", 0)

    monkeypatch.setattr(pubchempy, "get_compounds",
                        lambda q, ns: (_ for _ in ()).throw(
                            pubchempy.NotFoundError(404, "not found", [])))
    assert dock.resolve_cid("不存在的药") is None

    class BusyErr(pubchempy.PubChemHTTPError):
        def __init__(self):
            super().__init__(503, "Server Busy", [])

    monkeypatch.setattr(pubchempy, "get_compounds",
                        lambda q, ns: (_ for _ in ()).throw(BusyErr()))
    with pytest.raises(dock.PubChemUnavailable, match="503"):
        dock.resolve_cid("metformin")


def test_resolve_cid_success(monkeypatch):
    """正常解析:返回首个化合物的 CID。"""
    import pubchempy
    monkeypatch.setattr(pubchempy, "get_compounds",
                        lambda q, ns: [type("C", (), {"cid": 4091})()])
    assert dock.resolve_cid("metformin") == 4091


def test_download_pdb_caches_and_reuses(tmp_path, monkeypatch):
    """受体下载:首次联网(打桩)落盘;二次直接复用本地文件零请求。"""
    import requests

    calls = {"n": 0}

    def fake_get(url, timeout=None, **kw):
        calls["n"] += 1
        return _FakeResp("ATOM      1  CA  ALA\n")

    monkeypatch.setattr(requests, "get", fake_get)
    p1 = dock.download_pdb("2WXW")
    assert calls["n"] == 1 and open(p1, encoding="utf-8").read().startswith("ATOM")
    p2 = dock.download_pdb("2WXW")   # 缓存命中
    assert p2 == p1 and calls["n"] == 1


def test_smiles_to_sdf_3d_real_rdkit(tmp_path):
    """SMILES → 3D SDF:RDKit 纯本地;非法 SMILES 报错。"""
    out = dock.smiles_to_sdf_3d("CCO", str(tmp_path / "lig.sdf"), n_confs=3)
    assert out.endswith("lig.sdf") and open(out).read().count("END") > 0
    with pytest.raises(ValueError):
        dock.smiles_to_sdf_3d("not-a-smiles", str(tmp_path / "bad.sdf"))


def test_run_obabel_missing_tool_friendly_error(tmp_path, monkeypatch):
    """obabel 缺失:报错信息含工具路径指引,不是裸 FileNotFoundError。"""
    monkeypatch.setattr(config, "OBABEL_EXE", "")
    with pytest.raises(RuntimeError, match="[Oo]pen [Bb]abel|OBABEL"):
        dock._run_obabel(str(tmp_path / "in.sdf"), str(tmp_path / "out.pdbqt"))


# ---------- geo:E-utilities 检索 ----------


def test_search_geo_datasets_parses(monkeypatch):
    """GEO 检索:esearch 取 UID 列表 → esummary 解析标题/样本数/平台;细胞系关键词标注。"""
    import pandas as pd

    class Resp:
        def __init__(self, payload):
            self._p = payload

        def json(self):
            return self._p

    def fake_request(url, **kw):
        if "esearch" in url:
            return Resp({"esearchresult": {"idlist": ["200000001", "200000002"]}})
        return Resp({"result": {"uids": ["200000001", "200000002"], "200000001": {
            "accession": "GSE1297", "title": "Alzheimer hippocampus",
            "n_samples": 20, "gpl": "GPL570"}, "200000002": {
                "accession": "GSE000CELL", "title": "HepG2 cell line treatment",
                "n_samples": 6, "gpl": "GPL6947"}}})

    monkeypatch.setattr(http_client, "request_with_retry", fake_request)
    df = geo.search_geo_datasets("alzheimer", limit=5)
    assert list(df["gse"]) == ["GSE1297", "GSE000CELL"]
    cell_flags = df.set_index("gse")["possible_cellline"]
    assert cell_flags["GSE000CELL"] == "是" and cell_flags["GSE1297"] == ""
    links = df.set_index("gse")["link"]
    assert links["GSE1297"].endswith("acc.cgi?acc=GSE1297")


def test_search_geo_datasets_empty(monkeypatch):
    """检索零命中:返回带列的空表(不崩溃)。"""
    class Resp:
        def json(self):
            return {"esearchresult": {"idlist": []}}

    monkeypatch.setattr(http_client, "request_with_retry", lambda url, **kw: Resp())
    df = geo.search_geo_datasets("no_such_disease")
    assert isinstance(df, pd.DataFrame) and df.empty
