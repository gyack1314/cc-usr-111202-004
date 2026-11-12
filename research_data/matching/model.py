"""资料模型：严格校验、精确解析、候选边构建。

身份规则、分层字段、缺失口径、caliper、权重与容量均在此处统一处理，
求解器只接收已经建好的候选边，不重复解释原始资料。
"""
from dataclasses import dataclass, field
from fractions import Fraction

from ..io import integer
from .numbers import exact_number

# 未建边原因代码（诊断输出使用）
REASON_STRATUM = 'stratum_mismatch'
REASON_MISSING_CASE = 'missing_covariate_case'
REASON_MISSING_CONTROL = 'missing_covariate_control'
REASON_CALIPER = 'caliper_violation'


@dataclass(frozen=True)
class Subject:
    id: str
    stratum: str
    covariates: dict          # 名称 -> Fraction 或 None
    outcome: Fraction | None


@dataclass(frozen=True)
class Edge:
    case_id: str
    control_id: str
    stratum: str
    standardized: dict        # 名称 -> 有符号标准化差（Fraction）
    cost: Fraction            # 配置权重下的平方标准化距离之和
    distance_exact: dict      # 名称 -> 绝对标准化距离（Fraction）


@dataclass
class Material:
    cases: list
    controls: list
    scales: dict
    weights: dict
    calipers: dict            # caliper 集合 id -> {维度: Fraction}
    k: int
    dims: list = field(default_factory=list)
    case_index: dict = field(default_factory=dict)
    control_index: dict = field(default_factory=dict)
    capacities: dict = field(default_factory=dict)

    def control_capacity(self):
        return sum(self.capacities.values())


def _rational_param(x, name):
    if isinstance(x, bool) or not isinstance(x, (str, int)):
        raise ValueError(f'{name} 需为整数或分数字符串')
    try:
        value = Fraction(x)
    except (ValueError, ZeroDivisionError):
        raise ValueError(f'{name} 不是有效有理数') from None
    return value


def parse_material(v):
    """解析并整体校验一份资料字典，返回 Material。"""
    if not isinstance(v, dict):
        raise ValueError('顶层需要对象')
    if not isinstance(v.get('cases'), list) or not v['cases']:
        raise ValueError('cases 需要非空列表')
    if not isinstance(v.get('controls'), list) or not v['controls']:
        raise ValueError('controls 需要非空列表')

    def parse_subject(row, side):
        if not isinstance(row, dict):
            raise ValueError(f'{side} 记录需要对象')
        sid = row.get('id')
        if not isinstance(sid, str) or not sid.strip():
            raise ValueError(f'{side} 身份不能为空')
        stratum = row.get('stratum')
        if not isinstance(stratum, str) or not stratum.strip():
            raise ValueError(f'{sid} 分层缺失')
        cov = row.get('covariates')
        if not isinstance(cov, dict) or not cov:
            raise ValueError(f'{sid} 协变量缺失')
        outcome = row.get('outcome')
        if outcome is None:
            outcome_fraction = None
        else:
            outcome_fraction = exact_number(outcome)
        parsed = {}
        for name, value in cov.items():
            parsed[name] = None if value is None else exact_number(value)
        return Subject(sid, stratum, parsed, outcome_fraction)

    cases = [parse_subject(r, '病例') for r in v['cases']]
    controls = [parse_subject(r, '对照') for r in v['controls']]

    case_ids = {s.id for s in cases}
    control_ids = {s.id for s in controls}
    if len(case_ids) != len(cases):
        raise ValueError('病例身份重复')
    if len(control_ids) != len(controls):
        raise ValueError('对照身份重复')
    if case_ids & control_ids:
        raise ValueError('病例与对照身份不能重叠')

    scales_raw = v.get('scales')
    if not isinstance(scales_raw, dict) or not scales_raw:
        raise ValueError('scales 需要非空对象')
    dims = list(scales_raw.keys())
    scales = {}
    for name in dims:
        s = exact_number(scales_raw[name])
        if s <= 0:
            raise ValueError(f'scale {name} 必须为外部锁定的正数')
        scales[name] = s

    for s in cases + controls:
        if set(s.covariates) != set(dims):
            raise ValueError(f'{s.id} 协变量列与 scales 不一致')

    weights_raw = v.get('cost_weights')
    if not isinstance(weights_raw, dict) or set(weights_raw) != set(dims):
        raise ValueError('cost_weights 需覆盖全部协变量维度')
    weights = {}
    for name in dims:
        w = _rational_param(weights_raw[name], f'cost_weights.{name}')
        if w < 0:
            raise ValueError(f'cost_weights.{name} 不能为负')
        weights[name] = w
    if not any(w > 0 for w in weights.values()):
        raise ValueError('cost_weights 至少一个维度为正')

    caliper_sets_raw = v.get('caliper_sets')
    if not isinstance(caliper_sets_raw, list) or not caliper_sets_raw:
        raise ValueError('caliper_sets 需要非空列表')
    caliper_sets = {}
    for cs in caliper_sets_raw:
        if not isinstance(cs, dict) or not isinstance(cs.get('id'), str) or not cs['id'].strip():
            raise ValueError('caliper 集合需要字符串 id')
        cid = cs['id']
        if cid in caliper_sets:
            raise ValueError(f'caliper 集合 id 重复: {cid}')
        body = {k: cs[k] for k in cs if k != 'id'}
        if set(body) != set(dims):
            raise ValueError(f'caliper 集合 {cid} 需覆盖全部维度')
        calipers = {}
        for name in dims:
            bound = _rational_param(body[name], f'caliper.{cid}.{name}')
            if bound < 0:
                raise ValueError(f'caliper.{cid}.{name} 不能为负')
            calipers[name] = bound
        caliper_sets[cid] = calipers

    k = integer(v.get('k'), 1)
    capacities = {}
    for s in controls:
        # 从原始字典读取，保证整数口径校验
        raw = next(r for r in v['controls'] if r['id'] == s.id)
        capacities[s.id] = integer(raw['capacity'], 1)

    return Material(
        cases=cases,
        controls=controls,
        scales=scales,
        weights=weights,
        calipers=caliper_sets,
        k=k,
        dims=dims,
        case_index={s.id: s for s in cases},
        control_index={s.id: s for s in controls},
        capacities=capacities,
    )


def evaluate_pair(material, case, control, calipers):
    """返回 (edge 或 None, 原因代码 或 None)。

    分层不同、任一侧协变量缺失（不按零距离处理）、任一维超出 caliper
    都不能建边。
    """
    if case.stratum != control.stratum:
        return None, REASON_STRATUM
    standardized = {}
    distance = {}
    for name in material.dims:
        cv = case.covariates[name]
        tv = control.covariates[name]
        if cv is None:
            return None, REASON_MISSING_CASE
        if tv is None:
            return None, REASON_MISSING_CONTROL
        z = (cv - tv) / material.scales[name]
        standardized[name] = z
        d = abs(z)
        distance[name] = d
        if d > calipers[name]:
            return None, REASON_CALIPER
    cost = sum(material.weights[name] * standardized[name] ** 2
               for name in material.dims)
    edge = Edge(case.id, control.id, case.stratum,
                standardized, cost, distance)
    return edge, None


def build_edges(material, caliper_set_id):
    """在指定 caliper 集合下穷举所有合法候选边并给出排除统计。"""
    if caliper_set_id not in material.calipers:
        raise ValueError(f'未知 caliper 集合: {caliper_set_id}')
    calipers = material.calipers[caliper_set_id]
    edges = []
    exclusions = []
    counts = {REASON_STRATUM: 0, REASON_MISSING_CASE: 0,
              REASON_MISSING_CONTROL: 0, REASON_CALIPER: 0}
    for case in material.cases:
        for control in material.controls:
            edge, reason = evaluate_pair(material, case, control, calipers)
            if edge is not None:
                edges.append(edge)
            else:
                counts[reason] += 1
                detail = {'case_id': case.id, 'control_id': control.id, 'reason': reason}
                if reason == REASON_CALIPER:
                    detail['violations'] = [
                        name for name in material.dims
                        if case.covariates[name] is not None
                        and control.covariates[name] is not None
                        and abs((case.covariates[name] - control.covariates[name])
                                / material.scales[name]) > calipers[name]
                    ]
                exclusions.append(detail)
    edges.sort(key=lambda e: (e.case_id, e.control_id))
    return edges, {'counts': counts, 'details': exclusions}
