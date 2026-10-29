"""pytest 插件:封外网跑测试(审计工具,不进默认流程)。

用途:证明测试套件零外网依赖,或拦截新测试里的意外联网:
    pytest -p WY.tools.net_block

拦截非回环 TCP 连接(connect/connect_ex/create_connection)与非本机名 DNS
解析(getaddrinfo/gethostbyname),加载即对整个测试进程生效(含 job 后台线程)。

放行 127.0.0.0/8 与 ::1:asyncio 在 Windows 建事件循环时用回环 socketpair
做自唤醒管道(进程内 IPC,不是联网),TestClient 每次实例化都会触发,不能拦。

边界:拦的是 Python socket 层,不是安全沙箱——subprocess 联网、raw socket
等绕过途径不设防。定位是"无意联网立刻打红"的回归防线,2026-10 封网审核
(全量 344 条在该插件下全绿)的产物。
"""
import socket

_LOOPBACK_HOSTS = {"localhost", "127.0.0.1", "::1", "0.0.0.0", ""}

_real_connect = socket.socket.connect
_real_connect_ex = socket.socket.connect_ex
_real_create_connection = socket.create_connection
_real_getaddrinfo = socket.getaddrinfo
_real_gethostbyname = socket.gethostbyname


def _is_loopback(address) -> bool:
    host = str(address[0] if isinstance(address, tuple) else address).lower()
    return host in _LOOPBACK_HOSTS or host.startswith("127.")


def _deny_connect(self, address):
    if _is_loopback(address):
        return _real_connect(self, address)
    raise OSError(f"NETWORK BLOCKED: connect({address!r})")


def _deny_connect_ex(self, address):
    if _is_loopback(address):
        return _real_connect_ex(self, address)
    raise OSError(f"NETWORK BLOCKED: connect_ex({address!r})")


def _deny_create_connection(address, *args, **kwargs):
    if _is_loopback(address):
        return _real_create_connection(address, *args, **kwargs)
    raise OSError(f"NETWORK BLOCKED: create_connection({address!r})")


def _deny_getaddrinfo(host, *args, **kwargs):
    if _is_loopback(host):
        return _real_getaddrinfo(host, *args, **kwargs)
    raise OSError(f"NETWORK BLOCKED: dns resolve({host!r})")


def _deny_gethostbyname(host):
    if _is_loopback(host):
        return _real_gethostbyname(host)
    raise OSError(f"NETWORK BLOCKED: dns resolve({host!r})")


def pytest_configure(config):
    socket.socket.connect = _deny_connect
    socket.socket.connect_ex = _deny_connect_ex
    socket.create_connection = _deny_create_connection
    socket.getaddrinfo = _deny_getaddrinfo
    socket.gethostbyname = _deny_gethostbyname
