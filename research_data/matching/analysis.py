"""匹配结果的复核分析：诊断、复用集中度、效应汇总与协变量平衡。

口径承诺：
- 效应只使用拿满 k 个对照的完整病例；部分病例单列，不静默丢弃；
- 每个完整病例先对自己的对照结果取均值，再算病例减对照，最后病例等权；
- 复用只报告次数与集中度，不把匹配边当作独立样本给标准误；
- 平衡表只描述已观测协变量的加权差异，不宣称消除未测混杂。
"""
from fractions import Fraction

from .model import (REASON_CALIPER, REASON_MISSING_CASE,
                    REASON_MISSING_CONTROL, REASON_STRATUM)

REASON_TEXT = {
    REASON_STRATUM: '分层字段不同',
    REASON_MISSING_CASE: '病例侧协变量缺失（不按零距离处理）',
    REASON_MISSING_CONTROL: '对照侧协变量缺失（不按零距离处理）',
    REASON_CALIPER: '至少一维标准化距离超过 caliper',
}


def unmatched_diagnostics(material, edges, exclusion_info, solution):
    """按病例给出未匹配/未满配原因，并附最小割瓶颈线索。"""
    eligible = {}
    for e in edges:
        eligible.setdefault(e.case_id, []).append(e.control_id)
    chosen = {}
    for cid, tid in solution['pairs']:
        chosen.setdefault(cid, []).append(tid)

    reachable = set(solution['min_cut']['reachable_from_source'])
    per_case = []
    for case in material.cases:
        matched = sorted(chosen.get(case.id, []))
        candidates = sorted(eligible.get(case.id, []))
        status = 'complete' if len(matched) == material.k else (
            'partial' if matched else 'unmatched')
        reasons = []
        if not candidates:
            detail = _dominant_exclusion(exclusion_info, case.id)
            reasons.append({'code': detail[0], 'text': REASON_TEXT[detail[0]],
                            'detail': detail[1]})
        elif len(matched) < material.k:
            cap_left = sum(
                material.capacities[c]
                - sum(1 for x, t in solution['pairs'] if t == c)
                for c in candidates)
            reasons.append({
                'code': 'capacity_or_competition',
                'text': '存在合法候选边，但全局容量/竞争约束下槽位不足',
                'detail': {'eligible_controls': candidates,
                           'residual_capacity_among_eligible': cap_left,
                           'needed': material.k - len(matched)},
            })
        in_cut = f'case:{case.id}' in reachable
        per_case.append({
            'case_id': case.id,
            'stratum': case.stratum,
            'status': status,
            'matched_controls': matched,
            'matched_count': len(matched),
            'eligible_control_count': len(candidates),
            'eligible_controls': candidates,
            'on_source_side_of_mincut': in_cut,
            'reasons': reasons,
        })

    bottleneck_controls = sorted(
        name.split(':', 1)[1] for name in reachable
        if name.startswith('control:'))
    return {
        'per_case': per_case,
        'exclusion_totals': exclusion_info['counts'],
        'exclusion_reason_text': REASON_TEXT,
        'mincut_source_side': sorted(reachable),
        'bottleneck_controls_in_cut': bottleneck_controls,
        'note': '最小割源侧病例即被容量瓶颈阻挡的病例；割上的对照->汇点弧'
                '对应被用满的对照容量。',
    }


def _dominant_exclusion(exclusion_info, case_id):
    rows = [d for d in exclusion_info['details'] if d['case_id'] == case_id]
    priority = [REASON_MISSING_CASE, REASON_MISSING_CONTROL,
                REASON_STRATUM, REASON_CALIPER]
    for code in priority:
        hit = [d for d in rows if d['reason'] == code]
        if hit:
            return code, {'blocked_pairs': len(hit),
                          'examples': [d['control_id'] for d in hit[:5]],
                          **({'violations_example': hit[0].get('violations')}
                             if code == REASON_CALIPER else {})}
    return REASON_CALIPER, {'blocked_pairs': 0, 'examples': []}


def reuse_and_concentration(material, solution):
    """对照复用次数与匹配槽位集中度（HHI 与最大份额）。"""
    total_slots = max(1, len(solution['pairs']))
    rows = []
    for control in material.controls:
        used = solution['control_flows'][control.id]
        cases = sorted(cid for cid, tid in solution['pairs']
                       if tid == control.id)
        rows.append({
            'control_id': control.id,
            'capacity': material.capacities[control.id],
            'times_used': used,
            'distinct_cases': len(cases),
            'cases': cases,
            'slot_share': Fraction(used, total_slots),
            'capacity_utilization': Fraction(
                used, material.capacities[control.id]),
        })
    shares = [Fraction(r['times_used'], total_slots) for r in rows]
    hhi = sum(s * s for s in shares)
    return {
        'per_control': rows,
        'total_matched_slots': len(solution['pairs']),
        'reused_controls': [r['control_id'] for r in rows
                            if r['times_used'] > 1],
        'max_slot_share': max(shares),
        'herfindahl_hirschman_index': hhi,
        'warning': '复用边不是独立观测；本模块不提供把全部匹配边当独立'
                   '样本的标准误。',
    }


def effect_summary(material, solution):
    """完整病例等权的病例减对照差；部分病例单列。"""
    by_case = {}
    for cid, tid in solution['pairs']:
        by_case.setdefault(cid, []).append(tid)
    complete, partial = [], []
    diffs = []
    for case in material.cases:
        controls = sorted(by_case.get(case.id, []))
        record = {'case_id': case.id, 'matched_controls': controls,
                  'matched_count': len(controls)}
        if len(controls) == material.k:
            if case.outcome is None:
                raise ValueError(f'完整病例 {case.id} 缺少 outcome')
            control_outcomes = []
            for tid in controls:
                o = material.control_index[tid].outcome
                if o is None:
                    raise ValueError(f'对照 {tid} 缺少 outcome')
                control_outcomes.append(o)
            control_mean = sum(control_outcomes) / len(control_outcomes)
            diff = case.outcome - control_mean
            diffs.append(diff)
            record.update({'case_outcome': case.outcome,
                           'control_outcome_mean': control_mean,
                           'case_minus_controls': diff})
            complete.append(record)
        else:
            record['note'] = '未拿满 k 个对照，不进入效应汇总，仅单列'
            partial.append(record)
    pooled = sum(diffs) / len(diffs) if diffs else None
    return {
        'estimand': '完整病例等权：每病例对照结果均值 -> 病例减对照 -> 均值',
        'complete_cases': complete,
        'partial_or_unmatched_cases': partial,
        'n_complete_cases': len(complete),
        'n_excluded_partial_cases': len(partial),
        'effect_estimate': pooled,
        'standard_error': None,
        'standard_error_note': '按设计不提供把匹配边当独立样本的标准误；'
                               '推断需另用与复用结构相容的方法。',
    }


def _weighted_mean(subjects, weights, dim):
    active = [(s, weights[s.id]) for s in subjects
              if s.covariates[dim] is not None and weights[s.id] > 0]
    total_w = sum(w for _, w in active)
    if total_w == 0:
        return None, 0
    mean = sum(s.covariates[dim] * w for s, w in active) / total_w
    return mean, total_w


def balance_table(material, solution):
    """匹配前后按病例权重与对照使用权重计算的协变量差异（精确分数）。"""
    pre_case_w = {s.id: 1 for s in material.cases}
    pre_control_w = {s.id: 1 for s in material.controls}
    post_case_w = {s.id: 0 for s in material.cases}
    post_control_w = {s.id: 0 for s in material.controls}
    for cid, tid in solution['pairs']:
        post_case_w[cid] += 1
        post_control_w[tid] += 1

    rows = []
    for dim in material.dims:
        row = {'dimension': dim, 'scale': material.scales[dim]}
        for phase, cw, tw in (('pre_match', pre_case_w, pre_control_w),
                              ('post_match', post_case_w, post_control_w)):
            cm, cn = _weighted_mean(material.cases, cw, dim)
            tm, tn = _weighted_mean(material.controls, tw, dim)
            diff = None if cm is None or tm is None else cm - tm
            smd = None if diff is None else diff / material.scales[dim]
            row[phase] = {
                'case_weighted_mean': cm,
                'control_weighted_mean': tm,
                'case_weight_n': cn,
                'control_weight_n': tn,
                'mean_difference': diff,
                'standardized_mean_difference': smd,
            }
        rows.append(row)
    return {
        'dimensions': rows,
        'weights_rule': '匹配前各主体权重 1；匹配后病例权重=其匹配槽数，'
                        '对照权重=其被使用次数（容量内复用按次数计权）。',
        'missing_rule': '缺失协变量不参与该维度均值，权重人数单独报告，'
                        '缺失绝不按零距离进入匹配。',
        'interpretation': '平衡仅针对已观测协变量；不能解读为未测混杂已消除。',
    }
