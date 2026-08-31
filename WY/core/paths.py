"""项目路径锚点。

core 包内的模块不再位于 WY 根目录，原先依赖
`os.path.dirname(os.path.abspath(__file__))` 定位 Static 缓存目录的写法
在迁移后会把路径错指到 WY/core/ 下。统一改用 BASE_DIR（WY 包根目录）锚定，
保证拆分后磁盘上的缓存路径与旧版完全一致。
"""

import os

# WY 包根目录（core/ 的上一级）
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
