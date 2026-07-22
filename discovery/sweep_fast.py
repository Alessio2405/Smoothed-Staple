"""Fast rigorous-ish sweep: grid-optimized n=8 polyline per beta, then
exact-hull-margin scaling correction (reported L is achievable)."""
import numpy as np, json, sys
from fast import setup, optimize as gopt
from explore import margins as exact_margins

ZETA = 2.2782916414368529
def comps(beta):
    t = np.tan(beta)
    return dict(square=3*t*(2*np.cos(beta))/(t+2),
                zee=6*np.cos(beta)/np.sqrt(9+1/t**2),
                caliper=ZETA*np.sin(beta), diam=2*np.cos(beta))

P = np.array(json.load(open('polished.json'))['8']['P'])
P -= P.mean(0)
rho0 = np.sin(2*np.deg2rad(36.0))
angles = [float(a) for a in sys.argv[2:]]
rows=[]
for bdeg in angles:
    beta = np.deg2rad(bdeg)
    rho = np.sin(2*beta)
    seed = (P * (rho/rho0)).flatten()
    S = setup(beta, ndirs=2880)
    L, Pn, mmg = gopt(beta, 8, seed, S, maxiter=600)
    # exact correction: scale to restore exact feasibility
    mg, cyc = exact_margins(beta, Pn)
    # per-edge distances
    dmin = None
    # recompute distances from margins: D = cr^2 - rho^2 e^2 -> dist = sqrt(D/e^2 + rho^2)... easier direct:
    from explore import W_points
    from scipy.spatial import ConvexHull
    pts, r2 = W_points(beta, Pn)
    cyc2 = list(ConvexHull(pts).vertices); Q = pts[cyc2]
    dists=[]
    for k in range(len(Q)):
        p,q = Q[k], Q[(k+1)%len(Q)]; e=q-p
        dists.append((p[0]*q[1]-p[1]*q[0])/np.linalg.norm(e))
    dmin = min(dists)
    s = max(1.0, rho/dmin)
    Lc = L*s
    c = comps(beta)
    print(f"b={bdeg:6.3f}  L={Lc:.8f} (scale corr {s-1:+.1e})  "
          f"sq {c['square']:.7f}  zee {c['zee']:.7f}  cal {c['caliper']:.7f}", flush=True)
    rows.append(dict(b=bdeg, L=float(Lc), **{k:float(v) for k,v in c.items()}))
    P = Pn - Pn.mean(0); rho0 = rho
json.dump(rows, open(sys.argv[1],'w'), indent=1)
