"""确定性合成资料生成器（线性同余，无第三方依赖）。

仅用于压力测试与演示；所有受试者与结果均为虚构。种子固定时输出逐字节
稳定，方便复核人员重复实验。生成的资料字典与 data/ 下 JSON 同构。
"""


class Lcg:
    def __init__(self, seed=20261112):
        self.state = seed & 0x7FFFFFFF

    def integer(self, low, high):
        # Numerical Recipes 常数
        self.state = (1664525 * self.state + 1013904223) & 0x7FFFFFFF
        return low + self.state % (high - low + 1)


def synthesize(n_cases=40, n_controls=30, n_dims=3, seed=20261112,
               missing_rate=0.08, k=3, strata=3):
    """生成同构资料；对照容量在 1..k 间，部分对照允许多病例复用。

    对照在随机选定的层内、围绕该层某个病例的协变量小扰动生成，使同一层
    内存在大量合法候选边与真实的容量竞争；部分主体随机缺失协变量。
    """
    rng = Lcg(seed)
    dims = [f'x{i + 1}' for i in range(n_dims)]
    scales = {d: 10 ** rng.integer(0, 1) for d in dims}
    weights = {d: str(rng.integer(1, 4)) for d in dims}
    stratum_of = [f's{rng.integer(0, strata - 1)}' for _ in range(n_cases)]

    def maybe_missing(value):
        return None if rng.integer(1, 100) <= int(missing_rate * 100) else value

    cases = []
    for i in range(n_cases):
        cases.append({
            'id': f'case{i}',
            'stratum': stratum_of[i],
            'covariates': {d: maybe_missing(rng.integer(20, 90)) for d in dims},
            'outcome': rng.integer(0, 1),
        })

    controls = []
    for i in range(n_controls):
        anchor = cases[rng.integer(0, n_cases - 1)]
        cov = {}
        for d in dims:
            base = anchor['covariates'][d]
            if base is None or rng.integer(1, 100) <= int(missing_rate * 100):
                cov[d] = None
            else:
                # 小扰动：标准化距离大多在 caliper 内，偶尔越界
                cov[d] = max(1, base + rng.integer(-12, 12))
        controls.append({
            'id': f'ctrl{i}',
            'stratum': anchor['stratum'],
            'covariates': cov,
            'capacity': rng.integer(1, k),
            'outcome': rng.integer(0, 1),
        })

    caliper_sets = [
        {'id': 'tight', **{d: '1/2' for d in dims}},
        {'id': 'wide', **{d: '3' for d in dims}},
    ]
    return {
        'cases': cases,
        'controls': controls,
        'scales': scales,
        'cost_weights': weights,
        'caliper_sets': caliper_sets,
        'k': k,
    }
