# 多中心病例的容量约束可比队列匹配

这是研究资料维护与全局队列匹配分析的独立Python项目。`data/`保存合成记录，`research_data.io`负责严格JSON读取与SHA-256来源摘要，`research_data.domain`提供当前原始资料统计，`research_data.cli`将摘要输出为JSON。字段、时间、计数和量纲见[资料约定](docs/data-contract.md)。所有研究、中心及受试者均为虚构。

`research_data.matching`子包提供容量约束下的**全局**队列匹配（非逐病例贪心）：三层词典序目标（最大槽数 -> 最小精确费用 -> 匹配边列表字典序）、最小割与费用最优性证书、独立核对入口、未匹配诊断、对照复用集中度、完整病例等权效应汇总、匹配前后平衡表与caliper扫描。方法学、精确费用口径与证书细节见[匹配方法学](docs/matching-methodology.md)。

Python 3.11及以上版本可直接使用，当前没有第三方依赖。在项目根目录运行：

```sh
python -I -m unittest discover -s tests -t .
python -I -m compileall -q research_data tests run_matching.py
python -I -m research_data.cli inspect data/normal.json
```

匹配分析通过隔离入口运行（`-I` 使解释器忽略数据目录，防止资料中的同名模块被导入；见方法学文档第7节）：

```sh
python -I run_matching.py match data/normal.json tight
python -I run_matching.py verify data/normal.json tight
python -I run_matching.py sweep data/normal.json
python -I run_matching.py synthesize --cases 60 --controls 45 --out /tmp/synthetic.json
```

`match`默认附带独立核对（结构、容量、最大流、最小割、精确费用、势能证书、双引擎重解、小资料穷举）。`--method`可选`staged`（Fraction精确引擎，默认）、`staged-int`（整数缩放引擎）与诊断性的`zero-cycle`；大资料穷举由`--exhaustive-limit`控制并自动跳过。

`inspect`只显示已有摘要，不修改来源。文件名以`invalid-`起始的资料包含故意不合法的输入，其他文件是同结构的业务边界。数据摘要与来源指纹应同时保留，原始JSON不覆盖。

可使用固定镜像容器：`docker build -t research-material .`，随后`docker run --rm research-material`查看资料；`docker run --rm research-material python -I -m unittest discover -s tests -t .`运行现有检查。全部资料在项目内，运行不需要外部账号或服务。
