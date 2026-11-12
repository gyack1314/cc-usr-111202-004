"""按预设 caliper 集合逐一重求全局解并比较。

每个 caliper 集合独立重建候选边与网络：缺失协变量在任何放宽口径下都
不能建边，更不会变成零距离。比较流量、费用与入选人群的增减。
"""
from .model import build_edges
from .solve import solve
from .certify import verify_solution


def run_sweep(material, method='auto', verify=True):
    results = []
    previous = None
    for cid in material.calipers:
        edges, exclusions = build_edges(material, cid)
        solution = solve(material, edges, method=method)
        entry = {
            'caliper_set': cid,
            'calipers': material.calipers[cid],
            'candidate_edges': len(edges),
            'exclusion_totals': exclusions['counts'],
            'max_flow': solution['max_flow'],
            'fraction_total_cost': solution['fraction_total_cost'],
            'scaled_total_cost': solution['scaled_total_cost'],
            'pairs': solution['pairs'],
            'complete_cases': sorted(
                c for c, d in solution['case_flows'].items()
                if d == material.k),
            'used_controls': sorted(
                c for c, d in solution['control_flows'].items() if d > 0),
        }
        if verify:
            entry['verification'] = verify_solution(
                material, edges, solution, exhaustive_limit=64)['ok']
        if previous is not None:
            entry['change_vs_previous'] = _diff(previous, entry)
        results.append(entry)
        previous = entry
    return {
        'caliper_sets': results,
        'note': 'caliper 放宽只改变距离门槛；缺失协变量的排除在所有集合下'
                '一致，不会变成零距离。',
    }


def _diff(old, new):
    old_pairs = {tuple(p) for p in old['pairs']}
    new_pairs = {tuple(p) for p in new['pairs']}
    return {
        'flow_delta': new['max_flow'] - old['max_flow'],
        'candidate_edge_delta': new['candidate_edges'] - old['candidate_edges'],
        'pairs_added': sorted(new_pairs - old_pairs),
        'pairs_removed': sorted(old_pairs - new_pairs),
        'complete_cases_added': sorted(
            set(new['complete_cases']) - set(old['complete_cases'])),
        'complete_cases_lost': sorted(
            set(old['complete_cases']) - set(new['complete_cases'])),
        'controls_newly_used': sorted(
            set(new['used_controls']) - set(old['used_controls'])),
        'controls_released': sorted(
            set(old['used_controls']) - set(new['used_controls'])),
    }
