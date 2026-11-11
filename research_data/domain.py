"""读取病例协变量、外部标准化尺度及对照容量。"""
from .io import unique,positive,finite,integer

def summarize(v):
    cases=unique(v['cases'],'id');controls=unique(v['controls'],'id')
    if cases&controls:raise ValueError('双方身份不能重叠')
    scales=v['scales']
    if not scales:raise ValueError('缺少协变量尺度')
    for x in scales.values():positive(x)
    missing=0
    for row in v['cases']+v['controls']:
        if set(row['covariates'])!=set(scales):raise ValueError('协变量列不一致')
        if not isinstance(row['stratum'],str) or not row['stratum']:raise ValueError('分层缺失')
        for value in row['covariates'].values():
            if value is None:missing+=1
            else:finite(value)
    capacity=sum(integer(x['capacity'],1) for x in v['controls'])
    return {'cases':len(cases),'controls':len(controls),'control_capacity':capacity,'missing_cells':missing}
