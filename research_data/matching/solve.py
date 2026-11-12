"""三层目标全局求解器。

目标词典序（禁止用任意大罚数一次性凑合）：
1. 最大化总匹配槽数（最大流，允许病例只拿到 0..k 个对照）；
2. 在最大流解中最小化总费用（有理数精确费用的最小费用最大流）；
3. 在仍为最大流且费用最优的解中，最小化排序后的
   ``(case_id, control_id)`` 二元组列表的字典序。

实现方式（主路径，perturbation）：先求最大流 F；再把字典序编码成一组
扰动加到费用上，对固定流量 F 求一次最小费用流。扰动有显式整数上界
证明，保证它只能在“费用相同”的解之间决胜，绝不改变第 2 层的费用最优：

- 候选边按字典序编号 r=0..E-1，扰动 eps_r = -B^(E-1-r)，B=F+1；
  负的超递增权重使排名越靠前的边被优先选入（最小化扰动和）；
- 任一 F 条边的扰动绝对值之和 <= F·B^(E-1) < B^E；
- 整数缩放费用整体左移 E 位（乘 B^E），因此最小费用差（缩放后至少 1
  个整数单位）对应的 B^E 严格大于任何扰动差，费用层永远支配字典序层；
- Fraction 引擎把整数目标整体除以 lcm（正数保序），保持精确有理数，
  与整数引擎优化完全相同的全序。

这不是“调一个足够大的罚数”：B、E、F 全部由资料确定，支配关系可机器
核验。负权弧由 SPFA/Bellman-Ford 处理（初始网络为 DAG，无负环）。

提供两个互不共享中间量的数值引擎：
- ``staged``：Fraction 精确有理数 + SPFA；
- ``staged-int``：分母 LCM 等价整数缩放后的整数 + SPFA。
另有 ``zero-cycle`` 独立路径，在原始费用的最小费用流残量网络上做纯零
费用环强制，用于交叉核对第 3 层（其结果须与主路径一致）。
"""
from fractions import Fraction

from .network import (FlowGraph, augment_along, bfs_max_flow,
                      min_cost_flow_spfa, shortest_path_potential,
                      spfa_shortest_path)
from .numbers import common_denominator

S = 0


def node_layout(material):
    nc = len(material.cases)
    nt = len(material.controls)
    case_node = {s.id: 1 + i for i, s in enumerate(material.cases)}
    control_node = {s.id: 1 + nc + i for i, s in enumerate(material.controls)}
    t = 1 + nc + nt
    return case_node, control_node, t


def build_flow_graph(material, edges, cost_of):
    """构建 S->病例->对照->T 网络；cost_of(edge) 返回弧费用。"""
    case_node, control_node, t = node_layout(material)
    g = FlowGraph(t + 1)
    for s in material.cases:
        g.add_edge(S, case_node[s.id], material.k, 0)
    for s in material.controls:
        g.add_edge(control_node[s.id], t, material.capacities[s.id], 0)
    meta = []
    by_pair = {}
    for e in edges:
        arc = g.add_edge(case_node[e.case_id], control_node[e.control_id],
                         1, cost_of(e), ('edge', e.case_id, e.control_id))
        meta.append(arc)
        by_pair[(e.case_id, e.control_id)] = arc
    return g, meta, by_pair, t


def maximum_flow(material, edges):
    """第一层：Edmonds-Karp 最大流与最小割证书。"""
    g, meta, _, t = build_flow_graph(material, edges, lambda e: 0)
    flow, reachable = bfs_max_flow(g, S, t)
    cut_arcs = []
    cut_capacity = 0
    for e in g.arcs():
        if e.rev is None:
            continue
        is_original = e.meta is not None or e.u == S or e.v == t
        if not is_original:
            continue
        if e.u in reachable and e.v not in reachable:
            original_cap = e.cap + e.rev.cap
            cut_arcs.append({
                'u': _node_name(material, e.u),
                'v': _node_name(material, e.v),
                'capacity': original_cap,
                'kind': _arc_kind(e),
            })
            cut_capacity += original_cap
    return flow, {
        'reachable_from_source': sorted(_node_name(material, v)
                                        for v in reachable),
        'cut_arcs': cut_arcs,
        'cut_capacity': cut_capacity,
    }


def _arc_kind(arc):
    if arc.meta is not None:
        return 'candidate_edge'
    if arc.u == S:
        return 'source_case'
    return 'control_sink'


def _node_name(material, node):
    case_node, control_node, t = node_layout(material)
    if node == S:
        return 'SOURCE'
    if node == t:
        return 'SINK'
    for cid, idx in case_node.items():
        if idx == node:
            return f'case:{cid}'
    for cid, idx in control_node.items():
        if idx == node:
            return f'control:{cid}'
    raise AssertionError('未知节点')


def integer_scaling(edges):
    """所有候选边 Fraction 费用的公共分母 LCM，费用乘 L 后恰为整数。"""
    lcm = common_denominator([e.cost for e in edges])
    scaled = {(e.case_id, e.control_id): e.cost * lcm for e in edges}
    for value in scaled.values():
        assert value.denominator == 1
    return lcm, {k: int(v) for k, v in scaled.items()}


def perturbation_params(edges, required):
    """返回 (B, E, B^E, 排名, 每条边的负扰动, 扰动绝对值和上界)。

    负超递增扰动 eps_r = -B^(E-1-r)：最小化扰动和等价于最小化排序后
    边列表的字典序。支配界与穷举可核验。
    """
    e_count = len(edges)
    base = required + 1
    cost_place = base ** e_count
    rank = {(e.case_id, e.control_id): r for r, e in enumerate(edges)}
    epsilon = {pair: -(base ** (e_count - 1 - r))
               for pair, r in rank.items()}
    # 任一可行解至多含 required 条边，扰动绝对值和上界 < B^E
    max_perturbation = required * (base ** max(e_count - 1, 0)) \
        if e_count else 0
    assert e_count == 0 or max_perturbation < cost_place
    return base, e_count, cost_place, rank, epsilon, max_perturbation


def solve_perturbed_fraction(material, edges, required, scaled_cost, lcm):
    """Fraction 引擎：整数扰动目标整体除以 lcm，保持精确有理数。

    与整数引擎仅相差一个正的常数因子 1/lcm，因此全序完全相同。
    """
    base, e_count, cost_place, rank, epsilon, bound = \
        perturbation_params(edges, required)

    def adjusted(edge):
        pair = (edge.case_id, edge.control_id)
        integer_objective = scaled_cost[pair] * cost_place + epsilon[pair]
        return Fraction(integer_objective, lcm)

    g, meta, _, t = build_flow_graph(material, edges, adjusted)
    _, objective = min_cost_flow_spfa(g, S, t, required)
    pairs, _ = _extract_pairs(meta)
    evidence = {
        'engine': 'fraction_spfa',
        'base_B': base, 'edge_count_E': e_count,
        'cost_place_B_pow_E': cost_place,
        'max_perturbation_abs_sum_bound': bound,
        'dominance': 'max_perturbation_abs_sum_bound < cost_place_B_pow_E',
        'epsilon_rule': 'eps_r = -B^(E-1-r)，负超递增，最小化即字典序最小',
        'fraction_rule': '整数目标整体除以 lcm，与整数引擎同全序',
        'adjusted_objective': objective,
    }
    return pairs, evidence


def solve_perturbed_integer(material, edges, required, scaled_cost):
    """整数引擎：LCM 缩放费用左移 + 负整数扰动，一次整数最小费用流。"""
    base, e_count, cost_place, rank, epsilon, bound = \
        perturbation_params(edges, required)

    def adjusted(edge):
        pair = (edge.case_id, edge.control_id)
        return scaled_cost[pair] * cost_place + epsilon[pair]

    g, meta, _, t = build_flow_graph(material, edges, adjusted)
    _, objective = min_cost_flow_spfa(g, S, t, required)
    pairs, _ = _extract_pairs(meta)
    evidence = {
        'engine': 'integer_spfa',
        'base_B': base, 'edge_count_E': e_count,
        'cost_place_B_pow_E': cost_place,
        'max_perturbation_abs_sum_bound': bound,
        'dominance': 'max_perturbation_abs_sum_bound < cost_place_B_pow_E',
        'epsilon_rule': 'eps_r = -B^(E-1-r)，负超递增，最小化即字典序最小',
        'adjusted_objective_scaled': objective,
    }
    return pairs, evidence


def solve_zero_cycle_fraction(material, edges, required):
    """诊断性路径：原始费用最小费用流 + 贪心纯零费用环逐边强制。

    对当前不在流中的候选边 e=(i,j)：
    - 在不含源/汇、且避开已锁定弧的残量网络上找 j->i 最短路 P；
    - 若 P 的费用恰为 -cost(e)，则 P 与 e（i->j）构成零费用纯环：
      沿 P 增广 1 单位再把 e 置为在流中，总流量与总费用都不变；
    - 否则 e 被锁出。

    注意：这是贪心过程，早期锁出的边可能阻断后续字典序更优的环，因此
    不保证全局列表字典序最优。权威的第 3 层实现是带严格界证明的扰动法
    （solve_perturbed_*）。本函数仅用于小资料诊断与教学对照。
    """
    g, meta, _, t = build_flow_graph(material, edges, lambda e: e.cost)
    _, mc_cost = min_cost_flow_spfa(g, S, t, required)
    decisions = {'locked_in': 0, 'locked_out': 0}
    blocked_forward = set()    # 锁出边的正向弧
    blocked_reverse = set()    # 锁入边的反向弧（保护已锁定流）
    locked_out = []
    for arc in meta:
        if arc.rev.cap == 1:
            # 当前在流中：锁入，禁止以后经其反向弧把它换出
            blocked_reverse.add(id(arc.rev))
            decisions['locked_in'] += 1
            continue
        # e 不在流中（正向容量 1）。保持其正向弧开放作为环的返回段，
        # j->i 路径本身不会经过 i->j 的正向弧。
        blocked = blocked_forward | blocked_reverse
        dist, parent = spfa_shortest_path(
            g, arc.v, arc.u, blocked=blocked, exclude_nodes={S, t})
        if arc.u in dist and dist[arc.u] + arc.cost == 0:
            augment_along(parent, arc.u, 1)   # 沿 j->i 路径增广
            # 再把目标边 e 置为在流中（其容量恒为 1）
            assert arc.cap == 1 and arc.rev.cap == 0
            arc.cap = 0
            arc.rev.cap = 1
            blocked_reverse.add(id(arc.rev))
            decisions['locked_in'] += 1
        else:
            blocked_forward.add(id(arc))
            locked_out.append((arc.meta[1], arc.meta[2]))
            decisions['locked_out'] += 1
    pairs, _ = _extract_pairs(meta)
    assert len(pairs) == required
    return pairs, {'engine': 'fraction_zero_cycle_forcing',
                   'decisions': decisions,
                   'locked_out_edges': sorted(locked_out)}


def _extract_pairs(meta):
    pairs = []
    flows = {}
    for arc in meta:
        if arc.rev.cap == 1:
            pair = (arc.meta[1], arc.meta[2])
            pairs.append(pair)
            flows[pair] = 1
    pairs.sort()
    return pairs, flows


def _fraction_cost(pairs, edge_cost):
    total = Fraction(0)
    for pair in pairs:
        total += edge_cost[pair]
    return total


def _residual_potentials(g, t):
    """多源最短路对偶势能：残量正弧 reduced cost 全部 >= 0（Fraction）。"""
    return shortest_path_potential(g)


def solve(material, edges, method='auto'):
    """完整三层求解，返回结构化结果（含割与证书素材）。"""
    edge_cost = {(e.case_id, e.control_id): e.cost for e in edges}
    flow, cut = maximum_flow(material, edges)
    lcm, scaled_cost = integer_scaling(edges)

    if method == 'auto':
        method = 'staged'

    if flow == 0:
        return _empty_solution(material, edges, flow, cut, lcm, method)

    if method == 'staged':
        pairs, lex_evidence = solve_perturbed_fraction(
            material, edges, flow, scaled_cost, lcm)
    elif method == 'staged-int':
        pairs, lex_evidence = solve_perturbed_integer(
            material, edges, flow, scaled_cost)
    elif method == 'zero-cycle':
        pairs, lex_evidence = solve_zero_cycle_fraction(
            material, edges, flow)
    else:
        raise ValueError(f'未知求解方法: {method}')

    assert len(pairs) == flow
    total_cost = _fraction_cost(pairs, edge_cost)

    # 费用证书统一在完整候选网络上按最终流重建后独立计算：
    # 该流是全局最小费用流，当且仅当此残量网络无负费用环。
    used_potential = _fraction_potentials_for_flow(material, edges, pairs)
    case_flows, control_flows = _flow_counts(material, pairs)
    return {
        'method': method,
        'max_flow': flow,
        'pairs': pairs,
        'fraction_total_cost': total_cost,
        'scaled_total_cost': int(total_cost * lcm),
        'scaling': {'lcm': lcm,
                    'rule': '每条边整数费用 = Fraction费用 * lcm',
                    'exact': True},
        'case_flows': case_flows,
        'control_flows': control_flows,
        'min_cut': cut,
        'cost_certificate': _potential_certificate(material, used_potential),
        'lex_evidence': lex_evidence,
    }


def _fraction_potentials_for_flow(material, edges, pairs):
    """在完整 Fraction 网络上复刻给定流，求多源最短路对偶势能。"""
    g, meta, by_pair, t = build_flow_graph(material, edges, lambda e: e.cost)
    chosen = set(pairs)
    for arc in meta:
        pair = (arc.meta[1], arc.meta[2])
        if pair in chosen:
            assert arc.cap == 1 and arc.rev.cap == 0
            arc.cap = 0
            arc.rev.cap = 1
    return _residual_potentials(g, t)


def _flow_counts(material, pairs):
    case_flows = {s.id: 0 for s in material.cases}
    control_flows = {s.id: 0 for s in material.controls}
    for cid, tid in pairs:
        case_flows[cid] += 1
        control_flows[tid] += 1
    return case_flows, control_flows


def _potential_certificate(material, potentials):
    case_node, control_node, t = node_layout(material)
    names = ['SOURCE'] + [f'case:{s.id}' for s in material.cases] \
        + [f'control:{s.id}' for s in material.controls] + ['SINK']
    return {
        'potentials': {names[v]: potentials[v] for v in range(t + 1)},
        'rule': '残量正弧 reduced_cost = cost + p[u] - p[v] 必须全部 >= 0'
                '（多源最短路对偶可行势能）',
    }


def _empty_solution(material, edges, flow, cut, lcm, method):
    zero = Fraction(0)
    return {
        'method': method,
        'max_flow': 0,
        'pairs': [],
        'fraction_total_cost': zero,
        'scaled_total_cost': 0,
        'scaling': {'lcm': lcm,
                    'rule': '每条边整数费用 = Fraction费用 * lcm',
                    'exact': True},
        'case_flows': {s.id: 0 for s in material.cases},
        'control_flows': {s.id: 0 for s in material.controls},
        'min_cut': cut,
        'cost_certificate': _potential_certificate(
            material, {v: zero for v in range(len(material.cases)
                                              + len(material.controls) + 2)}),
        'lex_evidence': {'engine': 'empty',
                         'decisions': {'locked_in': 0, 'locked_out': 0},
                         'locked_out_edges': []},
    }
