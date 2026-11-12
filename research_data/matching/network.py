"""流网络与基础算法（全部标准库实现，支持 Fraction 或 int 费用）。

提供三套互相独立的底层算法，供求解与独立核对分别使用：
- ``bfs_max_flow``：Edmonds-Karp 广度优先最大流（只看容量不看费用）；
- ``spfa_shortest_path``：Bellman-Ford/SPFA 最短路，接受负边与 Fraction；
- ``dijkstra_shortest_path``：势能 + Dijkstra，要求输入势能使简化费用非负。

残量网络约定：正向弧索引为偶数，反向弧索引为奇数，互为 ``rev``。
"""
import heapq
from collections import deque


class Edge:
    __slots__ = ('u', 'v', 'cap', 'cost', 'rev', 'meta')

    def __init__(self, u, v, cap, cost, meta=None):
        self.u = u
        self.v = v
        self.cap = cap
        self.cost = cost
        self.rev = None
        self.meta = meta          # ('edge', case_id, control_id) 或 None


class FlowGraph:
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]

    def add_edge(self, u, v, cap, cost, meta=None):
        fwd = Edge(u, v, cap, cost, meta)
        rev = Edge(v, u, 0, -cost, None)
        fwd.rev = rev
        rev.rev = fwd
        self.g[u].append(fwd)
        self.g[v].append(rev)
        return fwd

    def arcs(self):
        for adj in self.g:
            for e in adj:
                yield e


def bfs_max_flow(graph, s, t):
    """Edmonds-Karp 最大流，返回 (流量, 残余网络上 s 可达点集合)。"""
    flow = 0
    while True:
        parent = {s: None}
        queue = deque([s])
        while queue and t not in parent:
            u = queue.popleft()
            for e in graph.g[u]:
                if e.cap > 0 and e.v not in parent:
                    parent[e.v] = e
                    queue.append(e.v)
        if t not in parent:
            break
        # 找瓶颈
        aug = None
        e = parent[t]
        while e is not None:
            aug = e.cap if aug is None else min(aug, e.cap)
            e = parent[e.u]
        e = parent[t]
        while e is not None:
            e.cap -= aug
            e.rev.cap += aug
            e = parent[e.u]
        flow += aug
    reachable = set()
    queue = deque([s])
    reachable.add(s)
    while queue:
        u = queue.popleft()
        for e in graph.g[u]:
            if e.cap > 0 and e.v not in reachable:
                reachable.add(e.v)
                queue.append(e.v)
    return flow, reachable


def spfa_shortest_path(graph, s, t, blocked=None, exclude_nodes=None):
    """SPFA 求 s->t 最短路（允许负简化费用）。

    blocked: 被锁定排除的正向弧集合，这些弧视为容量 0。
    exclude_nodes: 这些节点禁止进入路径（零费用环搜索时排除源/汇，
    保证找到的是流守恒的纯环而非改变总流量的增广路）。
    返回 (距离字典, 前驱弧字典)；不可达时距离中无 t。
    """
    blocked = blocked or frozenset()
    exclude = frozenset(exclude_nodes or ())
    dist = {s: 0}
    parent = {s: None}
    in_queue = {s}
    queue = deque([s])
    while queue:
        u = queue.popleft()
        in_queue.discard(u)
        for e in graph.g[u]:
            if e.cap <= 0 or id(e) in blocked or e.v in exclude:
                continue
            nd = dist[u] + e.cost
            if e.v not in dist or nd < dist[e.v]:
                dist[e.v] = nd
                parent[e.v] = e
                if e.v not in in_queue:
                    in_queue.add(e.v)
                    queue.append(e.v)
    return dist, parent


def shortest_path_potential(graph):
    """多源最短路势能：所有节点距离初始化为 0 后 Bellman-Ford 松弛。

    等价于引入一个到每个节点都有 0 费用弧的超级源，因此每个节点都有
    有限距离（通常 <= 0）。最小费用流的残量网络无负费用环，松弛收敛后
    对每条残量正弧满足 d[v] <= d[u] + cost，即
    cost + d[u] - d[v] >= 0，构成对偶可行势能。
    若 n-1 轮后仍可松弛，说明存在负费用环（流非最小费用）。
    """
    dist = [0] * graph.n
    for _ in range(graph.n - 1):
        changed = False
        for u in range(graph.n):
            for e in graph.g[u]:
                if e.cap <= 0:
                    continue
                nd = dist[u] + e.cost
                if nd < dist[e.v]:
                    dist[e.v] = nd
                    changed = True
        if not changed:
            break
    for u in range(graph.n):
        for e in graph.g[u]:
            if e.cap > 0 and dist[u] + e.cost < dist[e.v]:
                raise AssertionError('残量网络存在负费用环，流不是最小费用')
    return {v: dist[v] for v in range(graph.n)}


def dijkstra_shortest_path(graph, s, t, potential, blocked=None,
                           exclude_nodes=None):
    """势能 Dijkstra。要求对所有残量正弧 c(u,v)+p[u]-p[v] >= 0。

    exclude_nodes: 禁止进入路径的节点（零费用环搜索时排除源/汇）。
    返回 (简化距离字典, 前驱弧字典, 实际距离字典)。
    """
    blocked = blocked or frozenset()
    exclude = frozenset(exclude_nodes or ())
    inf = None
    dist = {s: 0}
    parent = {s: None}
    heap = [(0, s)]
    while heap:
        d, u = heapq.heappop(heap)
        if d != dist[u]:
            continue
        for e in graph.g[u]:
            if e.cap <= 0 or id(e) in blocked or e.v in exclude:
                continue
            reduced = e.cost + potential[u] - potential[e.v]
            if reduced < 0:
                raise AssertionError('势能不满足非负简化费用，Dijkstra 前提被破坏')
            nd = d + reduced
            if e.v not in dist or nd < dist[e.v]:
                dist[e.v] = nd
                parent[e.v] = e
                heapq.heappush(heap, (nd, e.v))
    actual = {v: dist[v] + potential[v] - potential[s] for v in dist}
    return dist, parent, actual


def min_cost_flow_spfa(graph, s, t, required, blocked=None):
    """反复 SPFA 增广，直至流量达到 required；达不到则抛出异常。"""
    flow = 0
    cost = 0
    while flow < required:
        dist, parent = spfa_shortest_path(graph, s, t, blocked)
        if t not in dist:
            raise AssertionError(f'无法达到所需流量 {required}，仅到 {flow}')
        aug = required - flow
        e = parent[t]
        while e is not None:
            aug = min(aug, e.cap)
            e = parent[e.u]
        e = parent[t]
        while e is not None:
            e.cap -= aug
            e.rev.cap += aug
            cost += e.cost * aug
            e = parent[e.u]
        flow += aug
    return flow, cost


def min_cost_flow_dijkstra(graph, s, t, required, potential, blocked=None):
    """反复势能 Dijkstra 增广并就地更新势能；返回 (流量, 实际费用, 势能)。"""
    flow = 0
    cost = 0
    while flow < required:
        dist, parent, _ = dijkstra_shortest_path(graph, s, t, potential, blocked)
        if t not in dist:
            raise AssertionError(f'无法达到所需流量 {required}，仅到 {flow}')
        for v, d in dist.items():
            potential[v] += d
        aug = required - flow
        e = parent[t]
        while e is not None:
            aug = min(aug, e.cap)
            e = parent[e.u]
        e = parent[t]
        while e is not None:
            e.cap -= aug
            e.rev.cap += aug
            cost += e.cost * aug
            e = parent[e.u]
        flow += aug
    return flow, cost, potential


def augment_along(parent, target, amount):
    """沿前驱弧字典增广 amount 单位流量。"""
    e = parent[target]
    while e is not None:
        e.cap -= amount
        e.rev.cap += amount
        e = parent[e.u]
