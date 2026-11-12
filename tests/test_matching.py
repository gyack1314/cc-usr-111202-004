"""全局队列匹配模块测试。

小资料穷举三层优先级（流量->费用->字典序）、容量瓶颈、等费用决胜、
缺失协变量口径与证书核对；大资料只验证全局解与证书，不枚举全部分配。
"""
import json
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from research_data.io import load
from research_data.matching.model import (parse_material, build_edges,
                                          evaluate_pair, REASON_MISSING_CASE,
                                          REASON_STRATUM)
from research_data.matching.solve import (solve, integer_scaling,
                                          perturbation_params,
                                          solve_perturbed_fraction,
                                          solve_perturbed_integer,
                                          maximum_flow)
from research_data.matching.certify import verify_solution, exhaustive_optimum
from research_data.matching.analysis import (effect_summary, balance_table,
                                             reuse_and_concentration,
                                             unmatched_diagnostics)
from research_data.matching.sweep import run_sweep
from research_data.matching.synthetic import synthesize
from research_data.matching.numbers import exact_number, common_denominator

DATA = Path(__file__).resolve().parents[1] / 'data'
METHODS = ['staged', 'staged-int', 'zero-cycle']


def load_material(name, caliper='tight'):
    value, _ = load(DATA / name)
    material = parse_material(value)
    edges, exclusions = build_edges(material, caliper)
    return material, edges, exclusions


class ExactNumberTests(unittest.TestCase):
    def test_decimal_text_roundtrip(self):
        self.assertEqual(exact_number(0.1), Fraction(1, 10))
        self.assertEqual(exact_number(2), Fraction(2))
        self.assertEqual(common_denominator(
            [Fraction(1, 4), Fraction(1, 6)]), 12)

    def test_reject_nonfinite_and_bool(self):
        import math
        for bad in (True, float('inf'), float('nan'), '0.1'):
            with self.assertRaises(ValueError):
                exact_number(bad)
        self.assertTrue(math.isfinite(0.0))


class ThreeLayerObjectiveTests(unittest.TestCase):
    """三份资料 x 两种 caliper x 三种方法，全部须与穷举一致。"""

    def test_all_materials_methods_match_exhaustive(self):
        for name in ('normal.json', 'capacity-bottleneck.json',
                     'equal-cost-competition.json'):
            for caliper in ('tight', 'wide'):
                material, edges, _ = load_material(name, caliper)
                brute = exhaustive_optimum(material, edges)
                for method in METHODS:
                    with self.subTest(material=name, caliper=caliper,
                                      method=method):
                        sol = solve(material, edges, method=method)
                        self.assertEqual(sol['max_flow'], brute['flow'])
                        self.assertEqual(sol['fraction_total_cost'],
                                         brute['cost'])
                        self.assertEqual([tuple(p) for p in sol['pairs']],
                                         brute['pairs'])
                        report = verify_solution(material, edges, sol)
                        self.assertTrue(report['ok'])

    def test_equal_cost_competition_lexicographic(self):
        # 等费用下应优先纳入字典序最小的 (x1,c1)
        material, edges, _ = load_material('equal-cost-competition.json')
        sol = solve(material, edges)
        pairs = [tuple(p) for p in sol['pairs']]
        self.assertIn(('x1', 'c1'), pairs)
        self.assertEqual(pairs, sorted(pairs))
        # 两个等费用竞争方案费用相同
        alt = {('x1', 'c2'), ('x2', 'c1'), ('x2', 'c2')}
        cost_chosen = sum(
            next(e.cost for e in edges
                 if (e.case_id, e.control_id) == p) for p in pairs)
        cost_alt = sum(
            next(e.cost for e in edges
                 if (e.case_id, e.control_id) == p) for p in alt)
        self.assertEqual(cost_chosen, cost_alt)
        self.assertLess(pairs, sorted(alt))

    def test_perturbation_dominance_bound(self):
        material, edges, _ = load_material('normal.json')
        flow, _ = maximum_flow(material, edges)
        base, e_count, cost_place, _, _, bound = perturbation_params(
            edges, flow)
        self.assertEqual(base, flow + 1)
        self.assertLess(bound, cost_place)

    def test_two_numeric_engines_agree(self):
        for name in ('normal.json', 'equal-cost-competition.json'):
            material, edges, _ = load_material(name)
            flow, _ = maximum_flow(material, edges)
            lcm, scaled = integer_scaling(edges)
            p_frac, _ = solve_perturbed_fraction(
                material, edges, flow, scaled, lcm)
            p_int, _ = solve_perturbed_integer(material, edges, flow, scaled)
            self.assertEqual(p_frac, p_int)
            # 权威解须与穷举三层最优一致
            self.assertEqual(p_frac, exhaustive_optimum(material, edges)['pairs'])


class CapacityBottleneckTests(unittest.TestCase):
    def test_partial_matching_is_global_not_greedy(self):
        material, edges, _ = load_material('capacity-bottleneck.json')
        sol = solve(material, edges)
        # 总槽数最大但没有任何病例拿满 k=2 以外的组合：x1 拿 2 个
        self.assertEqual(sol['max_flow'], 2)
        case_flows = sol['case_flows']
        self.assertEqual(case_flows['x1'], 2)
        self.assertEqual(case_flows['x2'], 0)
        # 容量上限被遵守
        for cid, used in sol['control_flows'].items():
            self.assertLessEqual(used, material.capacities[cid])
        # 最小割容量等于最大流
        self.assertEqual(sol['min_cut']['cut_capacity'], sol['max_flow'])

    def test_objective_maximizes_slots_not_complete_cases(self):
        # 目标第一层是总槽数（独立最大流），完整病例数不参与求解目标；
        # 效应分析才按“是否拿满 k”区分完整与部分病例。
        value = {
            'cases': [
                {'id': 'a', 'stratum': 's', 'covariates': {'x': 0},
                 'outcome': 1},
                {'id': 'b', 'stratum': 's', 'covariates': {'x': 0},
                 'outcome': 0},
            ],
            'controls': [
                {'id': 't1', 'stratum': 's', 'covariates': {'x': 0},
                 'capacity': 1, 'outcome': 0},
                {'id': 't2', 'stratum': 's', 'covariates': {'x': 0},
                 'capacity': 1, 'outcome': 1},
            ],
            'scales': {'x': 1},
            'cost_weights': {'x': '1'},
            'caliper_sets': [{'id': 'c', 'x': '1'}],
            'k': 2,
        }
        material = parse_material(value)
        edges, _ = build_edges(material, 'c')
        sol = solve(material, edges)
        # 第一层：槽数等于独立最大流（2），与穷举三层最优一致
        self.assertEqual(sol['max_flow'], exhaustive_optimum(material, edges)['flow'])
        flow, _ = maximum_flow(material, edges)
        self.assertEqual(sol['max_flow'], flow)
        # 部分病例被效应分析如实单列，而不是被求解器强行补成完整
        eff = effect_summary(material, sol)
        self.assertEqual(
            eff['n_complete_cases'] + eff['n_excluded_partial_cases'], 2)
        self.assertEqual(
            eff['n_complete_cases'],
            sum(1 for d in sol['case_flows'].values() if d == material.k))


class MissingCovariateTests(unittest.TestCase):
    def test_missing_is_excluded_not_zero_distance(self):
        value, _ = load(DATA / 'normal.json')
        material = parse_material(value)
        x3 = material.case_index['x3']
        c3 = material.control_index['c3']
        calipers = material.calipers['tight']
        edge, reason = evaluate_pair(material, x3, c3, calipers)
        self.assertIsNone(edge)
        self.assertEqual(reason, REASON_MISSING_CASE)
        # 任何 caliper 集合下缺失都不能建边
        for cid in material.calipers:
            edges, info = build_edges(material, cid)
            self.assertTrue(
                all(e.case_id != 'x3' for e in edges))
            self.assertGreaterEqual(
                info['counts'][REASON_MISSING_CASE], 1)

    def test_stratum_mismatch_excluded(self):
        value, _ = load(DATA / 'normal.json')
        material = parse_material(value)
        x1 = material.case_index['x1']
        c3 = material.control_index['c3']
        edge, reason = evaluate_pair(
            material, x1, c3, material.calipers['tight'])
        self.assertIsNone(edge)
        self.assertEqual(reason, REASON_STRATUM)


class IntegerScalingTests(unittest.TestCase):
    def test_scaled_costs_are_exact_integers(self):
        material, edges, _ = load_material('normal.json')
        lcm, scaled = integer_scaling(edges)
        for pair, value in scaled.items():
            self.assertIsInstance(value, int)
            original = next(e.cost for e in edges
                            if (e.case_id, e.control_id) == pair)
            self.assertEqual(Fraction(value, lcm), original)

    def test_fraction_and_integer_engine_costs_agree(self):
        material, edges, _ = load_material('equal-cost-competition.json')
        lcm, scaled = integer_scaling(edges)
        sol_frac = solve(material, edges, method='staged')
        sol_int = solve(material, edges, method='staged-int')
        self.assertEqual(sol_frac['pairs'], sol_int['pairs'])
        self.assertEqual(sol_int['scaled_total_cost'],
                         int(sol_frac['fraction_total_cost'] * lcm))


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.material, self.edges, self.excl = load_material('normal.json')
        self.sol = solve(self.material, self.edges)

    def test_effect_uses_only_complete_cases(self):
        eff = effect_summary(self.material, self.sol)
        # 只有 x1 拿满 k=2
        self.assertEqual(eff['n_complete_cases'], 1)
        self.assertEqual(
            [c['case_id'] for c in eff['complete_cases']], ['x1'])
        partial_ids = [c['case_id']
                       for c in eff['partial_or_unmatched_cases']]
        self.assertIn('x2', partial_ids)
        self.assertIn('x3', partial_ids)
        self.assertIsNone(eff['standard_error'])
        # 手工核算 x1 的效应
        x1 = self.material.case_index['x1']
        control_mean = Fraction(0 + 1, 2)   # c1=0, c2=1
        self.assertEqual(eff['effect_estimate'],
                         x1.outcome - control_mean)

    def test_reuse_concentration_reports_shares(self):
        reuse = reuse_and_concentration(self.material, self.sol)
        c2 = next(r for r in reuse['per_control'] if r['control_id'] == 'c2')
        self.assertEqual(c2['times_used'], 2)
        self.assertEqual(c2['distinct_cases'], 2)
        self.assertLessEqual(c2['times_used'],
                             self.material.capacities['c2'])
        self.assertIn('c2', reuse['reused_controls'])
        self.assertEqual(sum(r['times_used'] for r in reuse['per_control']),
                         self.sol['max_flow'])

    def test_balance_table_exact_and_labelled(self):
        bal = balance_table(self.material, self.sol)
        for row in bal['dimensions']:
            for phase in ('pre_match', 'post_match'):
                cell = row[phase]
                if cell['mean_difference'] is not None:
                    self.assertIsInstance(cell['mean_difference'], Fraction)
        self.assertIn('未测混杂', bal['interpretation'])

    def test_unmatched_diagnostics_explains_status(self):
        diag = unmatched_diagnostics(
            self.material, self.edges, self.excl, self.sol)
        by_id = {p['case_id']: p for p in diag['per_case']}
        self.assertEqual(by_id['x1']['status'], 'complete')
        self.assertEqual(by_id['x3']['status'], 'unmatched')
        self.assertTrue(by_id['x3']['reasons'])
        self.assertEqual(
            by_id['x3']['reasons'][0]['code'], REASON_MISSING_CASE)


class SweepTests(unittest.TestCase):
    def test_sweep_recomputes_global_solutions(self):
        material, _, _ = load_material('normal.json')
        report = run_sweep(material, verify=True)
        self.assertEqual(len(report['caliper_sets']),
                         len(material.calipers))
        for entry in report['caliper_sets']:
            self.assertTrue(entry['verification'])

    def test_widening_never_reduces_max_flow(self):
        value = synthesize(n_cases=30, n_controls=25, seed=7)
        material = parse_material(value)
        report = run_sweep(material, verify=False)
        flows = [e['max_flow'] for e in report['caliper_sets']]
        self.assertEqual(flows, sorted(flows))


class LargeMaterialTests(unittest.TestCase):
    """大资料：验证全局解与全部证书，但不枚举全部分配。"""

    def test_large_solve_and_certify_without_enumeration(self):
        # 确定性稠密资料：同层、同协变量 -> 完全二分图，边数稳定超穷举上限
        n_cases, n_controls = 20, 15
        value = {
            'cases': [
                {'id': f'a{i}', 'stratum': 's',
                 'covariates': {'x': 10, 'y': 20}, 'outcome': i % 2}
                for i in range(n_cases)],
            'controls': [
                {'id': f'b{i}', 'stratum': 's',
                 'covariates': {'x': 10, 'y': 20}, 'capacity': 1 + i % 3,
                 'outcome': i % 2}
                for i in range(n_controls)],
            'scales': {'x': 1, 'y': 1},
            'cost_weights': {'x': '1', 'y': '1'},
            'caliper_sets': [{'id': 'c', 'x': '1', 'y': '1'}],
            'k': 3,
        }
        material = parse_material(value)
        edges, _ = build_edges(material, 'c')
        self.assertEqual(len(edges), n_cases * n_controls)
        self.assertGreater(len(edges), 64)
        for method in ('staged', 'staged-int'):
            with self.subTest(method=method):
                sol = solve(material, edges, method=method)
                report = verify_solution(
                    material, edges, sol, exhaustive_limit=64)
                self.assertTrue(report['ok'])
                exhaustive = next(
                    c for c in report['checks']
                    if c['check'] == 'exhaustive_three_levels')
                self.assertIn('跳过', exhaustive['detail'])
        # 容量上界：总槽数不超过对照总容量，且受 k*病例数限制
        total_capacity = sum(material.capacities.values())
        self.assertLessEqual(sol['max_flow'],
                             min(total_capacity, n_cases * material.k))


class DeterminismTests(unittest.TestCase):
    def test_synthetic_is_seed_stable(self):
        a = synthesize(seed=123)
        b = synthesize(seed=123)
        self.assertEqual(
            json.dumps(a, sort_keys=True), json.dumps(b, sort_keys=True))

    def test_solve_is_deterministic(self):
        material, edges, _ = load_material('equal-cost-competition.json')
        first = solve(material, edges)['pairs']
        second = solve(material, edges)['pairs']
        self.assertEqual(first, second)


class IsolationTests(unittest.TestCase):
    """不可信数据目录中的同名模块不得被解释器加载。"""

    def test_planted_module_in_data_dir_not_imported(self):
        import subprocess
        import sys
        with tempfile.TemporaryDirectory() as d:
            dd = Path(d)
            (dd / 'study.json').write_text(
                (DATA / 'normal.json').read_text(encoding='utf8'),
                encoding='utf8')
            (dd / 'json.py').write_text(
                "raise RuntimeError('planted json module')\n",
                encoding='utf8')
            entry = Path(__file__).resolve().parents[1] / 'run_matching.py'
            proc = subprocess.run(
                [sys.executable, '-I', str(entry), 'match',
                 'study.json', 'tight', '--no-verify'],
                cwd=dd, capture_output=True, text=True, timeout=60)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            result = json.loads(proc.stdout)
            self.assertEqual(result['solution']['max_flow'], 3)


if __name__ == '__main__':
    unittest.main()
