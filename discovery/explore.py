"""Generalized escape-path search for T(beta), legs 1, base angle beta.

Escape criterion (Lutwak reduction, as in bellman-staple):
  path escapes  <=>  hull H satisfies  disk(rho) subset W(H),
  W(H) = Minkowski sum of w_i * R_{n_i} H over the 3 edge normals,
  rho = 2*Area = sin(2 beta).  Translation-invariant (sum w_i R_i = 0).

We minimize the length of a polyline through n points whose convex hull
must satisfy the escape condition. Path order fixed per run (permutation);
length = sum of consecutive distances. Hull = ConvexHull of all points.
"""
import itertools, json, sys
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

def triangle(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    nw = [(np.array([0.0,-1.0]), 2*cb),
          (np.array([ sb, cb ]), 1.0),
          (np.array([-sb, cb ]), 1.0)]
    return nw, 2*sb*cb

def rot(n):
    return np.array([[n[0], n[1]], [-n[1], n[0]]])

def W_points(beta, P):
    nw, rho = triangle(beta)
    copies = [w*(P @ rot(n).T) for n, w in nw]
    pts = np.array([a+b+c for a,b,c in itertools.product(*copies)])
    return pts, rho

def margins(beta, P, cycle=None):
    pts, rho = W_points(beta, P)
    # center W so that disk placement is optimal-ish: use Chebyshev-free version:
    # criterion is translation invariant in H, but W itself is fixed once H fixed
    # (origin fine since sum w_i R_i = 0 kills H-translation; W is centered by construction? not
    # necessarily, but the disk must sit at origin per the exact reduction used in the repo.)
    if cycle is None:
        cycle = list(ConvexHull(pts).vertices)
    Q = pts[cycle]
    n = len(Q)
    out = np.empty(n)
    for k in range(n):
        p, q = Q[k], Q[(k+1) % n]
        e = q - p
        cr = p[0]*q[1] - p[1]*q[0]
        out[k] = cr*cr - rho*rho*(e @ e)
        if cr < 0:  # wrong orientation / origin outside: strongly infeasible
            out[k] = -abs(out[k]) - 1.0
    return out, cycle

def path_length(P, order):
    Q = P[list(order)]
    return np.sum(np.linalg.norm(np.diff(Q, axis=0), axis=1))

def optimize_path(beta, n, order, x0, sym=False, iters=6, seed_jitter=0.0, rng=None):
    """P encoded as flat array. If sym: points come in mirror pairs/axis pts.
    Here sym=False: fully free 2n params."""
    x0 = np.asarray(x0, float).copy()
    if seed_jitter and rng is not None:
        x0 += rng.normal(0, seed_jitter, x0.shape)
    def unpack(x): return x.reshape(n, 2)
    _, cyc = margins(beta, unpack(x0))
    best = None
    for _ in range(iters):
        res = minimize(lambda x: path_length(unpack(x), order), x0,
                       constraints=[dict(type='ineq',
                                         fun=lambda x: margins(beta, unpack(x), cyc)[0])],
                       method='SLSQP', options=dict(ftol=1e-14, maxiter=600))
        x0 = res.x
        mg, cyc_new = margins(beta, unpack(x0))
        if cyc_new == cyc and mg.min() > -1e-11:
            best = (path_length(unpack(x0), order), unpack(x0).copy(), mg.min())
            break
        cyc = cyc_new
    else:
        mg, _ = margins(beta, unpack(x0))
        best = (path_length(unpack(x0), order), unpack(x0).copy(), mg.min())
    return best

def staple_seed():
    a, b, c = 0.2431, 0.18758, 0.45150
    # path: (-a,c) -> (-b,0) -> (b,0) -> (a,c)
    return np.array([[-a,c],[-b,0],[b,0],[a,c]])

if __name__ == '__main__':
    beta = np.deg2rad(36.0)
    rng = np.random.default_rng(0)

    # 0) sanity: reproduce the staple
    P0 = staple_seed()
    L, P, mm = optimize_path(beta, 4, (0,1,2,3), P0.flatten())
    print(f"[staple n=4 boundary] L = {L:.10f}  minmargin {mm:+.2e}")
    print("target (paper):        1.2849615334")

    results = [dict(tag='staple4', L=L, order=[0,1,2,3], P=P.tolist(), mm=mm)]

    # 1) n=4, all essentially-distinct Hamiltonian orders (zee-like uses diagonals)
    for order in [(1,0,2,3),(0,2,1,3),(1,0,3,2)]:
        bestL = np.inf
        for t in range(8):
            L, P, mm = optimize_path(beta, 4, order, P0.flatten(), seed_jitter=0.06*t/8+1e-9, rng=rng)
            if mm > -1e-9 and L < bestL:
                bestL, bP, bmm = L, P, mm
        print(f"[n=4 order {order}] L = {bestL:.10f}  mm {bmm:+.2e}")
        results.append(dict(tag=f'n4-{order}', L=bestL, order=list(order), P=bP.tolist(), mm=bmm))

    json.dump(results, open('explore_n4.json','w'), indent=1)
