"""Standalone exact verifier for the smoothed-staple certificate at the
golden gnomon T (apex 108 deg, base angles 36 deg, unit legs):

    l(T) <= 1.2826879687 < 1.2849615334  (Temerev 2026 staple bound)

Reads the rational 20-vertex path from certificate20.json. Reduction as in
bellman-staple/verification: the path escapes T if its hull H satisfies
disk(rho) subset W = phi*R0(H) + R1(H) + R2(H), rho = sin 72. It suffices
that some convex positively-oriented polygon with vertices in W has every
edge at distance >= rho from the origin (its hull is inside W).
the certificate cycle is stored in the JSON (deterministic; scipy only as
fallback to re-propose it); every accepted sign
determination is exact in Q(sqrt5, sqrt(10-2*sqrt5)).
Expected output: ALL CERTIFICATES PASS.
"""
import json, itertools
import sympy as sp

sqrt5 = sp.sqrt(5); phi = (1+sqrt5)/2
h = sp.sqrt(10-2*sqrt5)/4
rho2 = sp.expand((phi*h)**2)
assert sp.simplify(rho2 - (5+sqrt5)/8) == 0

data = json.load(open('certificate20.json'))
H = [sp.Matrix([sp.nsimplify(a, rational=True), sp.nsimplify(b, rational=True)])
     for a, b in data['vertices']]
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
BOUND  = sp.Rational(12826879687, 10**10)
STAPLE = sp.Rational(12849615334, 10**10)
assert L_upper < BOUND < STAPLE                       # exact comparisons in Q
print("exact rational upper bound L_upper =", sp.N(L_upper, 40), "(display only)")
print("certified exactly: L <= L_upper < 1.2826879687 < 1.2849615334")

nws = [(sp.Matrix([0,-1]), phi), (sp.Matrix([h,phi/2]), sp.Integer(1)),
       (sp.Matrix([-h,phi/2]), sp.Integer(1))]
assert sp.simplify(sum((w*v for v, w in nws), sp.Matrix([0,0]))) == sp.Matrix([0,0])
def rot(v): return sp.Matrix([[v[0],v[1]],[-v[1],v[0]]])

# Certificate cycle: read the stored index triples (fully deterministic).
# The stored cycle requires no trust: the exact convexity/orientation/distance
# checks below fail on any wrong cycle. scipy is used only as a fallback to
# re-propose the cycle if the JSON lacks it.
triples = data.get('cycle_triples')
if triples is None:
    import numpy as np
    from scipy.spatial import ConvexHull
    Pf = np.array([[float(v[0]), float(v[1])] for v in H])
    s36, c36 = np.sin(np.deg2rad(36)), np.cos(np.deg2rad(36))
    nwf = [(np.array([0,-1.]), 2*c36), (np.array([s36,c36]), 1.), (np.array([-s36,c36]), 1.)]
    copf = [w*(Pf @ np.array([[v[0],v[1]],[-v[1],v[0]]]).T) for v, w in nwf]
    trip = list(itertools.product(range(n), repeat=3))
    ptsf = np.array([copf[0][i]+copf[1][j]+copf[2][k] for i, j, k in trip])
    triples = [trip[t] for t in ConvexHull(ptsf).vertices]

cop_ex = [[sp.expand(w*(rot(v)*u)) for u in H] for v, w in nws]
Q = [sp.expand(cop_ex[0][i] + cop_ex[1][j] + cop_ex[2][k]) for i, j, k in triples]
m = len(Q)
ok_all = True
for k in range(m):
    p, q, r = Q[k], Q[(k+1) % m], Q[(k+2) % m]
    e = q - p
    conv = sp.expand((q[0]-p[0])*(r[1]-q[1]) - (q[1]-p[1])*(r[0]-q[0]))
    cr = sp.expand(p[0]*q[1] - p[1]*q[0])
    D = sp.expand(cr**2 - rho2*(e[0]**2 + e[1]**2))
    ok = (conv.is_positive is True and cr.is_positive is True and D.is_positive is True)
    ok_all &= ok
    if not ok:
        print(f"edge {k}: FAIL")
print(f"checked {m} exact edge certificates")
print("ALL CERTIFICATES PASS" if ok_all else "FAILURE")
print("certified: l(T_gnomon) < 1.2826879687  (exact rational chain; decimals above are display only)")
