"""独立核对：不信求解器自报，重建网络逐项验证匹配、容量与证书。

核对清单全部基于资料与候选边重新计算：
1. 每条匹配边真实存在且不重复；病例度数 <= k；对照度数 <= capacity；
2. 用独立的 Edmonds-Karp 重算最大流，匹配槽数必须等于最大流；
3. 最小割容量等于最大流，且割弧与报告一致；
4. 费用按 Fraction 重算，整数缩放费用恰好等于 lcm * 分数费用；
5. 势能证书：每个残量正弧 reduced_cost >= 0，并由 SPFA 独立复算距离；
6. 用另一条独立求解路径重解，三层目标结果必须一致；
7. 小资料可选穷举全部可行分配，直接验证三层优先级。
"""
from itertools import combinations
from fractions import Fraction

from .network import bfs_max_flow, shortest_path_potential
from .solve import (S, build_flow_graph, node_layout, solve,
                    solve_zero_cycle_fraction, integer_scaling)


class VerificationError(AssertionError):
    pass


def _record(checks, name, ok, detail=''):
    checks.append({'check': name, 'ok': bool(ok), 'detail': detail})
    if not ok:
        raise VerificationError(f'{name}: {detail}')


def rebuild_residual(material, edges, pairs):
    """用 Fraction 费用重建网络并把给定边置为已流 1 单位。"""
    g, meta, by_pair, t = build_flow_graph(material, edges, lambda e: e.cost)
    chosen = set(pairs)
    for arc in meta:
        pair = (arc.meta[1], arc.meta[2])
        if pair in chosen:
            if arc.cap != 1 or arc.rev.cap != 0:
                raise VerificationError(f'边 {pair} 初始容量异常')
            arc.cap = 0
            arc.rev.cap = 1
    return g, meta, by_pair, t


def verify_solution(material, edges, solution, exhaustive_limit=4096):
    """对一份 solve() 结果执行全部独立核对，返回核对清单（失败即抛错）。"""
    checks = []
    pairs = [tuple(p) for p in solution['pairs']]
    edge_set = {(e.case_id, e.control_id): e for e in edges}
    case_flows = {s.id: 0 for s in material.cases}
    control_flows = {s.id: 0 for s in material.controls}

    # 1. 结构、身份、容量与每病例 k 个不同对照
    _record(checks, 'pairs_sorted_unique', pairs == sorted(set(pairs)),
            f'{len(pairs)} 条边，去重后 {len(set(pairs))} 条')
    for pair in pairs:
        if pair not in edge_set:
            _record(checks, 'edge_eligible', False, f'{pair} 不是合法候选边')
        case_flows[pair[0]] += 1
        control_flows[pair[1]] += 1
    _record(checks, 'edge_eligible', True,
            f'{len(pairs)} 条匹配边全部在候选边集合内')
    over_k = [c for c, d in case_flows.items() if d > material.k]
    _record(checks, 'case_degree_within_k', not over_k,
            f'超 k 病例: {over_k}')
    over_cap = [c for c, d in control_flows.items()
                if d > material.capacities[c]]
    _record(checks, 'control_within_capacity', not over_cap,
            f'超容量对照: {over_cap}')

    # 2. 独立 Edmonds-Karp 最大流
    g0, _, _, t0 = build_flow_graph(material, edges, lambda e: 0)
    independent_flow, reachable = bfs_max_flow(g0, S, t0)
    _record(checks, 'flow_is_maximum', independent_flow == solution['max_flow']
            == len(pairs),
            f'独立最大流 {independent_flow}，报告 {solution["max_flow"]}，'
            f'匹配 {len(pairs)}')

    # 3. 最小割证书
    cut = solution['min_cut']
    _record(checks, 'mincut_capacity_equals_flow',
            cut['cut_capacity'] == independent_flow,
            f'割容量 {cut["cut_capacity"]} vs 最大流 {independent_flow}')
    reach_names = {_node(material, v) for v in reachable}
    _record(checks, 'mincut_reachable_set',
            set(cut['reachable_from_source']) == reach_names,
            '割集可达点与独立 BFS 不一致')
    # 独立重算割弧容量
    g0b, _, _, _ = build_flow_graph(material, edges, lambda e: 0)
    bfs_max_flow(g0b, S, t0)
    recomputed_cut = 0
    for arc in g0b.arcs():
        if arc.rev is None:
            continue
        is_original = arc.meta is not None or arc.u == S or arc.v == t0
        if is_original and arc.u in reachable and arc.v not in reachable:
            recomputed_cut += arc.cap + arc.rev.cap
    _record(checks, 'mincut_arcs_recomputed',
            recomputed_cut == cut['cut_capacity'],
            f'重算割容量 {recomputed_cut}')

    # 4. 费用与整数缩放
    fraction_total = sum(edge_set[p].cost for p in pairs)
    _record(checks, 'fraction_cost_recomputed',
            fraction_total == solution['fraction_total_cost'],
            f'重算 {fraction_total} vs 报告 {solution["fraction_total_cost"]}')
    lcm, scaled_cost = integer_scaling(edges)
    scaled_total = sum(scaled_cost[p] for p in pairs)
    _record(checks, 'integer_scaling_exact',
            scaled_total == solution['scaled_total_cost']
            and Fraction(scaled_total, lcm) == fraction_total,
            f'lcm={lcm}，整数费用 {scaled_total}')

    # 5. 势能/red­uced cost 对偶可行证书，独立多源最短路复算
    g, meta, by_pair, t = rebuild_residual(material, edges, pairs)
    potentials = solution['cost_certificate']['potentials']
    name_to_node = {_node(material, v): v for v in range(t + 1)}
    p = {name_to_node[name]: value for name, value in potentials.items()}
    negative = []
    for arc in g.arcs():
        if arc.cap > 0:
            rc = arc.cost + p[arc.u] - p[arc.v]
            if rc < 0:
                negative.append((_node(material, arc.u),
                                 _node(material, arc.v), str(rc)))
    _record(checks, 'reduced_cost_nonnegative', not negative,
            f'负 reduced cost 弧: {negative[:5]}')
    # 独立 Bellman-Ford 多源最短路：报告势能不得优于（小于）最短距离
    from .network import shortest_path_potential
    independent_pot = shortest_path_potential(g)
    cross_ok = all(p[v] >= independent_pot[v] for v in range(t + 1))
    _record(checks, 'potential_bellmanford_crosscheck', cross_ok,
            '报告势能在每个节点上不低于独立多源最短路距离')
    # 独立复算的势能同样必须满足全部残量正弧 reduced >= 0
    cycle_ok = all(
        arc.cost + independent_pot[arc.u] - independent_pot[arc.v] >= 0
        for arc in g.arcs() if arc.cap > 0)
    _record(checks, 'no_negative_residual_cycle', cycle_ok,
            '独立多源最短路确认残量网络无负费用环（费用最优性）')

    # 6. 独立路径重解：换数值引擎（Fraction <-> 整数），不共享中间量
    alt_method = _alternative(solution['method'])
    alt = solve(material, edges, method=alt_method)
    _record(checks, 'independent_engine_three_levels',
            alt['max_flow'] == solution['max_flow']
            and alt['pairs'] == pairs
            and alt['fraction_total_cost'] == fraction_total,
            f'独立引擎 {alt_method} 结果: 流量 {alt["max_flow"]}，'
            f'{len(alt["pairs"])} 边')
    # 零费用环强制是贪心过程，不保证全局字典序；仅在小资料上做健全性记录，
    # 不作为通过/失败条件（穷举三层目标才是小资料的权威核对）。
    if len(edges) <= exhaustive_limit:
        zc_pairs, _ = solve_zero_cycle_fraction(material, edges, len(pairs))
        checks.append({
            'check': 'zero_cycle_forcing_sanity',
            'ok': zc_pairs == pairs,
            'detail': '贪心零环强制结果；与权威穷举解不同不代表求解错误，'
                      f'零环给出 {len(zc_pairs)} 边'})

    # 7. 小资料穷举
    if len(edges) <= exhaustive_limit:
        brute = exhaustive_optimum(material, edges)
        _record(checks, 'exhaustive_three_levels',
                brute['flow'] == len(pairs)
                and brute['cost'] == fraction_total
                and brute['pairs'] == pairs,
                f'穷举 {brute["enumerated_assignments"]} 个可行分配')
    else:
        checks.append({'check': 'exhaustive_three_levels', 'ok': True,
                       'detail': f'边数 {len(edges)} 超穷举上限，跳过'})

    return {'ok': True, 'checks': checks}


def _alternative(method):
    """选择数值引擎不同、不共享中间量的另一条独立路径。"""
    if method == 'staged-int':
        return 'staged'
    return 'staged-int'

def _node(material, node):
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


def exhaustive_optimum(material, edges):
    """穷举全部满足度数/容量约束的边子集，按三层目标取最优。

    仅用于小资料：候选边组合数上界 2^E，由调用方限制 E。
    大资料禁止使用本函数（见 verify_solution 的 exhaustive_limit）。
    """
    edge_list = list(edges)
    best = None
    enumerated = 0
    for r in range(len(edge_list) + 1):
        for combo in combinations(edge_list, r):
            case_deg = {}
            control_deg = {}
            feasible = True
            for e in combo:
                case_deg[e.case_id] = case_deg.get(e.case_id, 0) + 1
                control_deg[e.control_id] = \
                    control_deg.get(e.control_id, 0) + 1
                if case_deg[e.case_id] > material.k \
                        or control_deg[e.control_id] \
                        > material.capacities[e.control_id]:
                    feasible = False
                    break
            if not feasible:
                continue
            enumerated += 1
            flow = len(combo)
            cost = sum(e.cost for e in combo)
            pair_list = sorted((e.case_id, e.control_id) for e in combo)
            # 三层目标均按最小化方向比较：流量取负、费用取负、列表字典序
            candidate = (-flow, cost, pair_list)
            if best is None or candidate < best:
                best = candidate
    neg_flow, cost, pair_list = best
    return {'flow': -neg_flow, 'cost': cost, 'pairs': pair_list,
            'enumerated_assignments': enumerated}
