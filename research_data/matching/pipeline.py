"""一次匹配分析的完整编排：解析 -> 建边 -> 求解 -> 证书 -> 复核分析。"""
from ..io import load
from .model import build_edges, parse_material
from .solve import solve
from .certify import verify_solution
from .analysis import (balance_table, effect_summary,
                       reuse_and_concentration, unmatched_diagnostics)


def run_file(path, caliper_set_id, method='auto', verify=True,
             exhaustive_limit=4096):
    value, digest = load(path)
    return run_material(value, digest, caliper_set_id, method, verify,
                        exhaustive_limit)


def run_material(value, digest, caliper_set_id, method='auto', verify=True,
                 exhaustive_limit=4096):
    material = parse_material(value)
    edges, exclusions = build_edges(material, caliper_set_id)
    solution = solve(material, edges, method=method)
    report = {
        'source_sha256': digest,
        'caliper_set': caliper_set_id,
        'calipers': material.calipers[caliper_set_id],
        'k': material.k,
        'dimensions': material.dims,
        'solution': solution,
        'candidate_edge_count': len(edges),
        'candidate_edges': [
            {'case_id': e.case_id, 'control_id': e.control_id,
             'cost': e.cost,
             'standardized': e.standardized,
             'distance_exact': e.distance_exact}
            for e in edges
        ],
        'unmatched_diagnostics': unmatched_diagnostics(
            material, edges, exclusions, solution),
        'reuse': reuse_and_concentration(material, solution),
        'effect': effect_summary(material, solution),
        'balance': balance_table(material, solution),
    }
    if verify:
        report['independent_verification'] = verify_solution(
            material, edges, solution, exhaustive_limit=exhaustive_limit)
    return material, edges, solution, report
