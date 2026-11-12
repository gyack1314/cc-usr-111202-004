"""命令行入口：match / verify / sweep / synthesize。

示例：
  python -I -m research_data.matching.cli match data/normal.json tight
  python -I -m research_data.matching.cli verify data/normal.json tight
  python -I -m research_data.matching.cli sweep data/normal.json
  python -I -m research_data.matching.cli synthesize --cases 60 --out /tmp/s.json
"""
import argparse
import sys

from ..io import load
from .model import build_edges, parse_material
from .numbers import dumps
from .pipeline import run_file
from .sweep import run_sweep
from .synthetic import synthesize
from .certify import VerificationError


def _write(text, out):
    if out:
        with open(out, 'w', encoding='utf8') as f:
            f.write(text)
    else:
        print(text)


def cmd_match(args):
    _, _, _, report = run_file(
        args.file, args.caliper_set, method=args.method,
        verify=not args.no_verify, exhaustive_limit=args.exhaustive_limit)
    _write(dumps(report), args.out)
    return 0


def cmd_verify(args):
    value, digest = load(args.file)
    material = parse_material(value)
    edges, _ = build_edges(material, args.caliper_set)
    from .solve import solve
    solution = solve(material, edges, method=args.method)
    from .certify import verify_solution
    result = verify_solution(material, edges, solution,
                             exhaustive_limit=args.exhaustive_limit)
    _write(dumps({'source_sha256': digest, 'verification': result}), args.out)
    return 0


def cmd_sweep(args):
    value, digest = load(args.file)
    material = parse_material(value)
    report = {'source_sha256': digest, 'sweep': run_sweep(
        material, method=args.method, verify=not args.no_verify)}
    _write(dumps(report), args.out)
    return 0


def cmd_synthesize(args):
    data = synthesize(n_cases=args.cases, n_controls=args.controls,
                      n_dims=args.dims, seed=args.seed, k=args.k)
    _write(dumps(data), args.out)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog='research_data.matching',
        description='容量约束全局队列匹配与独立证书核对')
    sub = parser.add_subparsers(dest='command', required=True)

    def add_common(p):
        p.add_argument('file')
        p.add_argument('caliper_set')
        p.add_argument('--method', choices=['auto', 'staged', 'staged-int', 'zero-cycle'],
                       default='auto')
        p.add_argument('--exhaustive-limit', type=int, default=4096)
        p.add_argument('--out')

    p_match = sub.add_parser('match', help='求全局匹配并输出完整分析')
    add_common(p_match)
    p_match.add_argument('--no-verify', action='store_true')
    p_match.set_defaults(func=cmd_match)

    p_verify = sub.add_parser('verify', help='独立核对入口')
    add_common(p_verify)
    p_verify.set_defaults(func=cmd_verify)

    p_sweep = sub.add_parser('sweep', help='逐 caliper 集合重求全局解')
    p_sweep.add_argument('file')
    p_sweep.add_argument('--method', choices=['auto', 'staged', 'staged-int', 'zero-cycle'],
                         default='auto')
    p_sweep.add_argument('--no-verify', action='store_true')
    p_sweep.add_argument('--out')
    p_sweep.set_defaults(func=cmd_sweep)

    p_syn = sub.add_parser('synthesize', help='生成确定性虚构资料')
    p_syn.add_argument('--cases', type=int, default=40)
    p_syn.add_argument('--controls', type=int, default=30)
    p_syn.add_argument('--dims', type=int, default=3)
    p_syn.add_argument('--k', type=int, default=3)
    p_syn.add_argument('--seed', type=int, default=20261112)
    p_syn.add_argument('--out')
    p_syn.set_defaults(func=cmd_synthesize)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (VerificationError, ValueError, KeyError, TypeError) as exc:
        print(f'核对失败或资料不合法: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
