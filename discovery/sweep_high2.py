import numpy as np, json, sys
sys.path.insert(0, 'bellman-staple/discovery')
from staple_family import polish as staple_polish
from fast import setup, optimize as gopt
from explore import W_points
from scipy.spatial import ConvexHull

ZETA = 2.2782916414368529
def comps(beta):
    t = np.tan(beta)
    return dict(square=3*t*(2*np.cos(beta))/(t+2),
                zee=6*np.cos(beta)/np.sqrt(9+1/t**2),
                caliper=ZETA*np.sin(beta))

def exact_correct(beta, Pn, L):
    pts, _ = W_points(beta, Pn)
    cyc = list(ConvexHull(pts).vertices); Q = pts[cyc]
    dmin = min((Q[k][0]*Q[(k+1)%len(Q)][1]-Q[k][1]*Q[(k+1)%len(Q)][0]) /
               np.linalg.norm(Q[(k+1)%len(Q)]-Q[k]) for k in range(len(Q)))
    s = max(1.0, np.sin(2*beta)/dmin)
    return L*s

def staple_seed_n8(beta):
    r = staple_polish(beta)
    a,b,c = r['a'], r['b'], r['c']
    V = np.array([[-a,c],[-b,0],[b,0],[a,c]])
    # subdivide legs into 3 pieces each -> 8 points
    def lerp(p,q,t): return p*(1-t)+q*t
    pts = [V[0], lerp(V[0],V[1],1/3), lerp(V[0],V[1],2/3), V[1], V[2],
           lerp(V[2],V[3],1/3), lerp(V[2],V[3],2/3), V[3]]
    return np.array(pts), r['L']

rows = []
prevP = None
for bdeg in [float(x) for x in sys.argv[1:]]:
    beta = np.deg2rad(bdeg)
    S = setup(beta, ndirs=2880)
    seeds = []
    Pst, Lst = staple_seed_n8(beta)
    seeds.append(Pst.flatten())
    if prevP is not None:
        seeds.append((prevP*np.sin(2*beta)/rho_prev).flatten())
    best = (np.inf, None)
    for sd in seeds:
        L, Pn, mm = gopt(beta, 8, sd, S, maxiter=600)
        Lc = exact_correct(beta, Pn, L)
        if Lc < best[0]: best = (Lc, Pn)
    Lc, Pn = best
    c = comps(beta)
    gain = Lst - Lc
    print(f"b={bdeg:6.3f}  L={Lc:.8f}  staple {Lst:.8f} (gain {gain:+.5f})  "
          f"sq {c['square']:.7f}  zee {c['zee']:.7f}", flush=True)
    rows.append(dict(b=bdeg, L=float(Lc), staple=float(Lst),
                     **{k:float(v) for k,v in c.items()}))
    prevP = Pn - Pn.mean(0); rho_prev = np.sin(2*beta)
json.dump(rows, open('sweep_high2.json','w'), indent=1)
