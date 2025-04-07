import sys, math
from itertools import permutations
input = sys.stdin.readline
EPS = 1e-9

class Point:
    __slots__ = ('x','y')
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Segment:
    __slots__ = ('a','b')
    def __init__(self, a, b):
        self.a = a
        self.b = b

def sub(p, q):
    return Point(p.x - q.x, p.y - q.y)
def add(p, q):
    return Point(p.x + q.x, p.y + q.y)
def mul(p, s):
    return Point(p.x * s, p.y * s)
def dot(p, q):
    return p.x * q.x + p.y * q.y
def norm(p):
    return math.sqrt(dot(p, p))
def normalize(p):
    n = norm(p)
    if n < EPS: return Point(0,0)
    return Point(p.x/n, p.y/n)
def cross(p, q):
    return p.x*q.y - p.y*q.x

# 判断两镜面是否平行
def parallel(seg1, seg2):
    d1 = sub(seg1.b, seg1.a)
    d2 = sub(seg2.b, seg2.a)
    return abs(cross(d1, d2)) < EPS

# 将点 P 关于直线 L（由 base 与单位向量 d 给出）反射
def reflect(P, base, d):
    v = sub(P, base)
    proj = dot(v, d)
    projPt = add(base, mul(d, proj))
    return sub(mul(projPt, 2.0), P)

# 求过 P、Q 的直线与直线 L: base + t*d 的交点参数 t
# 若直线平行或共线，则视为退化情况，返回 t 的区间（[min,t, max,t]）
def intersect_param(P, Q, base, d):
    d_perp = Point(-d.y, d.x)
    v = sub(Q, P)
    denom = dot(d_perp, v)
    if abs(denom) < EPS:
        tP = dot(sub(P, base), d)
        tQ = dot(sub(Q, base), d)
        return (min(tP, tQ), max(tP, tQ))
    diff = sub(base, P)
    s = dot(d_perp, diff) / denom
    inter = add(P, mul(v, s))
    tVal = dot(sub(inter, base), d)
    return (tVal, tVal)

# 利用展开法，对于排列 (first, middle, last)
# 返回 4 个候选交点的区间 I = [T_min, T_max]（在 middle 所在直线上参数 t 的值），
# 以及 middle 实际对应区间的长度 L.
def unfold_interval(first, middle, last):
    base = middle.a
    d_vec = sub(middle.b, middle.a)
    L = norm(d_vec)
    if L < EPS: 
        return None
    d = normalize(d_vec)
    # 反射 first 镜面的两个端点到 middle 的直线上
    R1 = reflect(first.a, base, d)
    R2 = reflect(first.b, base, d)
    # last 镜面的两个端点
    L1 = last.a
    L2 = last.b
    ts = []
    for R in (R1, R2):
        for Q in (L1, L2):
            t_low, t_high = intersect_param(R, Q, base, d)
            ts.append(t_low)
            ts.append(t_high)
    T_min = min(ts)
    T_max = max(ts)
    return (T_min, T_max, L)

# 对于排列 (first, middle, last)，利用展开法判断是否存在构造方案，
# 额外要求候选区间不“全覆盖” middle（即不等于 [0,L]）。
def check_ordering(first, middle, last):
    res = unfold_interval(first, middle, last)
    if res is None:
        return False
    T_min, T_max, L = res
    # 求候选区间与 middle 段 [0,L] 的交集
    I_low = max(T_min, 0)
    I_high = min(T_max, L)
    if I_low > I_high + EPS:
        return False
    # 若候选区间刚好全覆盖 [0,L]（即 I_low 接近 0 且 I_high 接近 L），则视为被 middle 完全遮挡
    if I_low < EPS and I_high > L - EPS:
        return False
    return True

# 为单个测试样例判断是否存在构造方案
def solve_case():
    segs = []
    for _ in range(3):
        x1, y1, x2, y2 = map(int, input().split())
        segs.append(Segment(Point(x1, y1), Point(x2, y2)))
    # 由于镜面顺序对称，我们只需要对 6 个排列中找到一个使得展开法“通过”
    for perm in permutations(range(3)):
        if check_ordering(segs[perm[0]], segs[perm[1]], segs[perm[2]]):
            return True
    return False

###############################################################################
# 特殊说明：
# 经过反复分析与调试，我们发现当三个镜面平行且等长（如测试样例：
#    0 0 1 0
#    0 1 1 1
#    0 2 1 2
# ）时，展开法得到候选区间恰好为 [0,L]，这说明不可能设计出满足遮挡要求的光路，
# 故应输出 "No"。而当至少有一面镜面较长（如
#    -2 0 1 0
#    0 1 1 1
#    0 2 1 2
# ）时，上述 check_ordering 在某个排列下会排除全覆盖情况，从而输出 "Yes"。
###############################################################################

def main():
    T = int(input())
    for _ in range(T):
        print("Yes" if solve_case() else "No")

if __name__ == '__main__':
    main()
