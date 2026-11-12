#!/usr/bin/env python3
"""隔离模式入口（请以 ``python -I run_matching.py ...`` 运行）。

数据文件路径始终通过命令行参数传入，解释器不会从数据所在目录加载模块，
因此资料目录中即使存在同名 json.py 等植入文件也不会被 import。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from research_data.matching.cli import main

if __name__ == '__main__':
    raise SystemExit(main())
