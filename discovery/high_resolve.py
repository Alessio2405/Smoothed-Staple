import numpy as np, json, sys
sys.path.insert(0,'bellman-staple/discovery')
from staple_family import polish as staple_polish
from fast import setup, optimize as gopt
from explore import W_points
from scipy.spatial import ConvexHull

def exact_corrected(beta, Pn, L):
    pts, _ = W_points(beta, Pn)
    Q = pts[list(ConvexHull(pts).vertices)]
    dmin = min((Q[k][0]*Q[(k+1)%len(Q)][1]-Q[k][1]*Q[(k+1)%len(Q)][0]) /
               np.linalg.norm(Q[(k+1)%len(Q)]-Q[k]) for k in range(len(Q)))
    return L*max(1.0, np.sin(2*beta)/dmin)

def staple_seed(beta, n=10):
    r = staple_polish(beta)
    a,b,c = r['a'], r['b'], r['c']
    V = np.array([[-a,c],[-b,0],[b,0],[a,c]])
    segs = [(V[0],V[1]),(V[1],V[2]),(V[2],V[3])]
    # distribute n points along the 3 segments proportional to length
    lens = [np.linalg.norm(q-p) for p,q in segs]
    ts = np.linspace(0, sum(lens), n)
    pts=[]; acc=0; k=0
    for t in ts:
        while k<2 and t>acc+lens[k]: acc+=lens[k]; k+=1
        f=(t-acc)/lens[k]; p,q=segs[k]
        pts.append(p*(1-f)+q*f)
    return np.array(pts), r['L']

prev=None
zee = lambda b: 6*np.cos(b)/np.sqrt(9+1/np.tan(b)**2)
out=[]
for bdeg in [42.0, 42.25, 42.5]:
    beta = np.deg2rad(bdeg)
    S = setup(beta, ndirs=2880)
    seeds=[]
    Pst, Lst = staple_seed(beta)
    seeds.append(Pst.flatten() + np.random.default_rng(1).normal(0,0.004,20))
    seeds.append(Pst.flatten())
    if prev is not None:
        seeds.append((prev*np.sin(2*beta)/rp).flatten())
    best=(np.inf,None)
    for sd in seeds:
        L,Pn,mm = gopt(beta, 10, sd, S, maxiter=800)
        Lc = exact_corrected(beta, Pn, L)
        if Lc < best[0]: best=(Lc,Pn)
    Lc,Pn = best
    z = zee(beta)
    print(f"b={bdeg:6.3f}: smoothed(n10,exact-corr) {Lc:.7f}  staple {Lst:.7f} "
          f"(gain {Lst-Lc:+.5f})  zee {z:.7f}  zee-margin {z-Lc:+.5f}", flush=True)
    out.append(dict(b=bdeg, L=float(Lc), staple=float(Lst), zee=float(z)))
    prev = Pn - Pn.mean(0); rp = np.sin(2*beta)

m = [(r['zee']-r['L']) for r in out]; bs=[r['b'] for r in out]
for i in range(len(bs)-1):
    if (m[i]>0)!=(m[i+1]>0):
        bc = bs[i]+(bs[i+1]-bs[i])*m[i]/(m[i]-m[i+1])
        print(f"\nsmoothed-vs-zee crossing (conservative): beta = {bc:.3f} deg  (staple: 42.28, Gibbs conj: 42.3)")
json.dump(out, open('high_resolve.json','w'), indent=1)
