"""Exact certificate from the 20-vertex path: cycle-only exact computation."""
import numpy as np, itertools, json
import sympy as sp
from scipy.spatial import ConvexHull

sqrt5 = sp.sqrt(5); phi = (1+sqrt5)/2
h = sp.sqrt(10-2*sqrt5)/4
rho2 = sp.expand((phi*h)**2)
assert sp.simplify(rho2 - (5+sqrt5)/8) == 0

P = np.load('P20.npy'); P -= P.mean(0)
DEN = 10**8; EPSF = 1e-6
Pr = [[sp.Rational(round(float(v)*(1+EPSF)*DEN), DEN) for v in row] for row in P]
H = [sp.Matrix(v) for v in Pr]
n = len(H)
L = sum(sp.sqrt((H[i+1][0]-H[i][0])**2+(H[i+1][1]-H[i][1])**2) for i in range(n-1))
Lv = sp.N(L, 40)
print("exact rational-path length L =", Lv)

nws = [(sp.Matrix([0,-1]), phi), (sp.Matrix([h,phi/2]), sp.Integer(1)),
       (sp.Matrix([-h,phi/2]), sp.Integer(1))]
def rot(nv): return sp.Matrix([[nv[0],nv[1]],[-nv[1],nv[0]]])

# float Minkowski with index bookkeeping -> propose cycle
Pf = np.array([[float(v[0]), float(v[1])] for v in Pr])
b36, c36f = np.sin(np.deg2rad(36)), np.cos(np.deg2rad(36))
nwf = [(np.array([0,-1.]), 2*c36f), (np.array([b36,c36f]),1.), (np.array([-b36,c36f]),1.)]
def rotf(nv): return np.array([[nv[0],nv[1]],[-nv[1],nv[0]]])
copf = [w*(Pf @ rotf(nv).T) for nv,w in nwf]
trip = list(itertools.product(range(n), repeat=3))
ptsf = np.array([copf[0][i]+copf[1][j]+copf[2][k] for i,j,k in trip])
cyc = list(ConvexHull(ptsf).vertices)
print("cycle size:", len(cyc))

# exact points only for the cycle
cop_ex = [[sp.expand(w*(rot(nv)*v)) for v in H] for nv,w in nws]
Q = []
for idx in cyc:
    i,j,k = trip[idx]
    Q.append(sp.expand(cop_ex[0][i]+cop_ex[1][j]+cop_ex[2][k]))
m = len(Q)
ok_all = True; worst = None
for k in range(m):
    p,q,r = Q[k], Q[(k+1)%m], Q[(k+2)%m]
    e = q-p
    conv = sp.expand((q[0]-p[0])*(r[1]-q[1])-(q[1]-p[1])*(r[0]-q[0]))
    cr = sp.expand(p[0]*q[1]-p[1]*q[0])
    D = sp.expand(cr**2 - rho2*(e[0]**2+e[1]**2))
    ok = (conv.is_positive is True and cr.is_positive is True and D.is_positive is True)
    ok_all &= ok
    d = float(sp.N(cr/sp.sqrt(e.dot(e)), 25))
    if worst is None or d < worst: worst = d
    if not ok: print(f" edge {k}: FAIL conv={conv.is_positive} cr={cr.is_positive} D={D.is_positive}")
print(f"min edge distance {worst:.10f} vs rho {float(sp.N(sp.sqrt(rho2),25)):.10f}"
      f"  margin {worst-float(sp.N(sp.sqrt(rho2),25)):+.2e}")
print("ALL CERTIFICATES PASS" if ok_all else "FAILURE")
print("\ncertified: l(T_gnomon) <=", sp.N(L, 30))
json.dump(dict(vertices=[[str(a),str(b)] for a,b in Pr], length_40d=str(Lv),
               cycle_size=m, pass_=bool(ok_all)), open('certificate20.json','w'), indent=1)
