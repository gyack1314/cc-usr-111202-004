"""精确有理数口径与 JSON 序列化辅助。

资料契约要求按 JSON 原十进制文本解析，不先转二进制浮点再声称有理精确。
Python 的 json 把数字解析为 int/float，对 float 取最短往返十进制文本
``str(x)`` 可还原写入时的十进制字面量（如 0.1 -> '0.1'），再经
``Decimal`` 转成精确分数。
"""
import json
import math
from decimal import Decimal
from fractions import Fraction


def exact_number(x):
    """有限十进制数值（int/float，拒绝 bool 与非有限值）转 Fraction。"""
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError('需要有限十进制数值')
    return Fraction(Decimal(str(x)))


def decimal_text(value, digits=10):
    """用 Decimal 做定点截断展示，绝不使用二进制浮点；返回字符串。"""
    if not isinstance(value, Fraction):
        raise TypeError('decimal_text 只接受 Fraction')
    quantum = Decimal(1).scaleb(-digits)
    text = str((Decimal(value.numerator) / Decimal(value.denominator)).quantize(quantum))
    return text


def common_denominator(fractions):
    """一组 Fraction 分母的最小公倍数；空序列返回 1。"""
    if not fractions:
        return 1
    d = 1
    for f in fractions:
        d = math.lcm(d, f.denominator)
    return d


class RationalEncoder(json.JSONEncoder):
    """Fraction 输出为分数字符串，tuple 输出为数组，禁用 NaN/Infinity。"""

    def default(self, o):
        if isinstance(o, Fraction):
            return str(o)
        if isinstance(o, tuple):
            return list(o)
        return super().default(o)


def dumps(obj, indent=2):
    return json.dumps(obj, ensure_ascii=False, indent=indent,
                      allow_nan=False, cls=RationalEncoder)
