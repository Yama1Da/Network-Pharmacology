"""自检脚本公共 runner:通过/失败统计与退出码,供 selfcheck_*.py 共用。"""

import os
import sys
import tempfile
from contextlib import contextmanager


def check(name: str, fn) -> bool:
    """运行单个测试并打印结果;断言失败记 FAIL,其他异常记 ERROR。"""
    try:
        fn()
        print(f"PASS  {name}")
        return True
    except AssertionError as e:
        print(f"FAIL  {name}: {e}")
        return False
    except Exception as e:
        print(f"ERROR {name}: {type(e).__name__}: {e}")
        return False


def run(tests: list[tuple[str, object]], title: str) -> None:
    """依次运行全部测试,全部通过时退出码为 0,有失败则为 1。"""
    print(f"==== {title}: {len(tests)} tests ====")
    ok = all([check(name, fn) for name, fn in tests])
    print("==== ALL PASS ====" if ok else "==== FAILED ====")
    sys.exit(0 if ok else 1)


@contextmanager
def patch_attr(obj, name: str, value):
    """临时替换对象属性,退出时恢复原值(测试打桩用)。"""
    orig = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield obj
    finally:
        setattr(obj, name, orig)


@contextmanager
def _tmp_uniprot_cache(config, name: str):
    """临时换装 UniProt 查询缓存:备份 config.UNIPROT_LOOKUP_CACHE → 清掉同名临时文件
    → 改指向临时文件;退出时恢复原值并删除临时文件。查询函数的替换请配合 patch_attr。"""
    tmp = os.path.join(tempfile.gettempdir(), name)
    if os.path.exists(tmp):
        os.remove(tmp)
    orig = config.UNIPROT_LOOKUP_CACHE
    config.UNIPROT_LOOKUP_CACHE = tmp
    try:
        yield tmp
    finally:
        config.UNIPROT_LOOKUP_CACHE = orig
        if os.path.exists(tmp):
            os.remove(tmp)
