"""自检脚本公共 runner:通过/失败统计与退出码,供 selfcheck_*.py 共用。"""

import sys


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
