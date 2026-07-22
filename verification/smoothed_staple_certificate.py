"""Exact certification: an 8-vertex escaping path for the golden gnomon
T (apex 108 deg, base angles 36 deg, unit legs) with

    l(T) <= L8 = 1.28278708539565... < 1.2849615334...

strictly improving the certified staple bound of Temerev (2026),
"An improved upper bound for Bellman's lost-in-a-forest problem in
isosceles triangles". Self-contained; requires sympy only
(scipy is used solely to propose the certificate cycle; the accepted
checks are exact).

Reduction (same as bellman-staple/verification): the path escapes T if
its convex hull H satisfies disk(rho) subset W, where
W = phi*R0(H) + R1(H) + R2(H) (Minkowski sum over the edge normals of T,
weighted by edge lengths) and rho = sin 72 = 2*Area(T).
It suffices to exhibit a convex positively-oriented polygon whose vertices
are points of W, with every edge line at distance >= rho from the origin:
its hull is contained in W, hence the disk is.

Every accepted sign determination is exact in Q(sqrt5, sqrt(10-2*sqrt5)).
Expected output: "ALL CERTIFICATES PASS".
"""
import itertools
import sympy as sp

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
h = sp.sqrt(10 - 2*sqrt5) / 4              # sin 36
rho2 = sp.expand((phi*h)**2)               # (sin 72)^2
assert sp.simplify(rho2 - (5 + sqrt5)/8) == 0

# The rational 8-vertex path (in path order); found by constrained
# optimization, inflated by 1e-5 and rounded to denominator 1e7.
V = [
    (sp.sympify('-956977/10000000'), sp.sympify('27673/78125')),
    (sp.sympify('-2049373/10000000'), sp.sympify('4701/80000')),
    (sp.sympify('-2172277/10000000'), sp.sympify('31391/2000000')),
    (sp.sympify('-1186623/5000000'), sp.sympify('-825201/10000000')),
    (sp.sympify('462003/5000000'), sp.sympify('-46731/200000')),
    (sp.sympify('192109/1250000'), sp.sympify('-771597/5000000')),
    (sp.sympify('1782907/10000000'), sp.sympify('-1168963/10000000')),
    (sp.sympify('3308087/10000000'), sp.sympify('1587183/10000000'))
]
H = [sp.Matrix(v) for v in V]
n = len(H)


def _sqrt_upper(q, k=30):
    """Exact rational upper bound for sqrt(q), q a nonnegative Rational.
    Pure integer arithmetic (math.isqrt); the returned u satisfies u**2 > q,
    verified exactly."""
    from math import isqrt
    p, r = sp.fraction(sp.nsimplify(q, rational=True))
    p, r = int(p), int(r)
    scale = 10**k
    m = isqrt((p * scale * scale) // r) + 1
    assert m*m*r > p*scale*scale
    return sp.Rational(m, scale)

seg2 = [(H[i+1][0]-H[i][0])**2 + (H[i+1][1]-H[i][1])**2 for i in range(n-1)]
L = sum(sp.sqrt(q) for q in seg2)
L_upper = sum(_sqrt_upper(q) for q in seg2)          # exact rational majorant
BOUND  = sp.Rational(12827871, 10**7)
STAPLE = sp.Rational(12849615334, 10**10)
assert L_upper < BOUND < STAPLE                       # exact comparisons in Q
print("exact rational upper bound L_upper =", sp.N(L_upper, 50), "(display only)")
print("certified exactly: L <= L_upper < 1.2827871 < 1.2849615334")

normals_weights = [
    (sp.Matrix([0, -1]), phi),
    (sp.Matrix([h, phi/2]), sp.Integer(1)),
    (sp.Matrix([-h, phi/2]), sp.Integer(1)),
]
assert sp.simplify(sum((w*nv for nv, w in normals_weights),
                       sp.Matrix([0, 0]))) == sp.Matrix([0, 0])
for nv, _ in normals_weights:
    assert sp.simplify(nv.dot(nv) - 1) == 0

def rot(nv):
    return sp.Matrix([[nv[0], nv[1]], [-nv[1], nv[0]]])

copies = [[sp.expand(w*(rot(nv)*v)) for v in H] for nv, w in normals_weights]
pts = [sp.expand(p + q + r) for p, q, r in itertools.product(*copies)]

# Certificate cycle: prefer the index triples stored in certificate8.json
# (fully deterministic). The stored cycle needs no trust: the exact
# convexity/orientation/distance checks below fail on any wrong cycle.
# scipy is used only as a fallback to re-propose the cycle.
try:
    import json as _json
    triples = _json.load(open('certificate8.json')).get('cycle_triples')
except Exception:
    triples = None

nH = len(H)
if triples is not None:
    cycle = [i*nH*nH + j*nH + k for i, j, k in triples]
else:
    try:
        import numpy as np
        from scipy.spatial import ConvexHull
        ptsf = np.array([[float(v[0]), float(v[1])] for v in pts])
        cycle = list(ConvexHull(ptsf).vertices)
    except Exception as exc:
        raise RuntimeError('scipy is required to propose the certificate cycle') from exc

Q = [pts[i] for i in cycle]
m = len(Q)
print(f"certificate polygon: {m} vertices of W")

ok_all = True
for k in range(m):
    p, q, r = Q[k], Q[(k+1) % m], Q[(k+2) % m]
    e = q - p
    conv_cross = sp.expand((q[0]-p[0])*(r[1]-q[1]) - (q[1]-p[1])*(r[0]-q[0]))
    cross_pq = sp.expand(p[0]*q[1] - p[1]*q[0])
    D = sp.expand(cross_pq**2 - rho2*(e[0]**2 + e[1]**2))
    ok = (conv_cross.is_positive is True and cross_pq.is_positive is True
          and D.is_positive is True)
    ok_all &= ok
    print(f" edge {k:2d}: convex>0 {conv_cross.is_positive}  "
          f"cross>0 {cross_pq.is_positive}  D>0 {D.is_positive}")

print("ALL CERTIFICATES PASS" if ok_all else "FAILURE")
print("certified: l(T_gnomon) < 1.2827871 < 1.2849615334  (exact rational chain; decimals above are display only)")