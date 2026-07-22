"""Generic piecewise line/arc path with template spec and constraint generation."""
import numpy as np
from scipy.optimize import minimize
np.seterr(all='ignore')

def tri(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    normals = np.array([[0,-1.],[sb,cb],[-sb,cb]])
    weights = np.array([2*cb,1.,1.])
    Rs = np.array([[[n[0],n[1]],[-n[1],n[0]]] for n in normals])
    return Rs, weights, 2*sb*cb

# template: list of piece kinds; params consumed in order
# 's'->1 param (length), 'a'->2 params (r, phi), 'c'->1 param (delta)
def build(template, x):
    pts=[np.zeros(2)]; arcs=[]; i=0
    p=np.zeros(2); th=x[0]; i=1
    for kind in template:
        if kind=='s':
            s=x[i]; i+=1
            p = p + s*np.array([np.cos(th),np.sin(th)]); pts.append(p.copy())
        elif kind=='a':
            r,ph = x[i],x[i+1]; i+=2
            if r>1e-9 and ph>1e-9:
                c = p + r*np.array([-np.sin(th),np.cos(th)])
                a0 = th - np.pi/2
                arcs.append((c.copy(), r, a0, a0+ph))
                p = c + r*np.array([np.cos(a0+ph),np.sin(a0+ph)])
                pts.append(p.copy())
            th += ph
        elif kind=='c':
            th += x[i]; i+=1
    return np.array(pts), arcs

def h_H(U, pts, arcs):
    h = (U @ pts.T).max(axis=1)
    g = np.arctan2(U[:,1], U[:,0])
    for c, r, a0, a1 in arcs:
        d = (g - a0) % (2*np.pi)
        inside = d <= (a1-a0)
        h = np.where(inside, np.maximum(h, U@c + r), h)
    return h

def margins(template, x, beta, U, cache):
    Rs, w, rho = cache
    pts, arcs = build(template, x)
    tot = np.zeros(len(U))
    for i in range(3):
        tot += w[i]*h_H(U @ Rs[i], pts, arcs)
    return tot - rho

def length(template, x):
    L=0.; i=1
    for kind in template:
        if kind=='s': L+=x[i]; i+=1
        elif kind=='a': L+=x[i]*x[i+1]; i+=2
        else: i+=1
    return L

def solve(template, x0, beta, rounds=14, m0=1024, verbose=False):
    cache = tri(beta)
    th = np.linspace(0,2*np.pi,m0,endpoint=False)
    U = np.c_[np.cos(th),np.sin(th)]
    thv = np.linspace(0,2*np.pi,400000,endpoint=False)
    Uv = np.c_[np.cos(thv),np.sin(thv)]
    # bounds
    bnds=[(None,None)]
    for kind in template:
        if kind=='s': bnds.append((0,2))
        elif kind=='a': bnds += [(0,6),(0,1.5)]
        else: bnds.append((0,2.6))
    x=np.asarray(x0,float)
    for it in range(rounds):
        res = minimize(lambda z: length(template,z), x,
                       constraints=[dict(type='ineq',
                                fun=lambda z: margins(template,z,beta,U,cache))],
                       bounds=bnds, method='SLSQP',
                       options=dict(ftol=1e-15, maxiter=2500))
        x = res.x
        mg = margins(template, x, beta, Uv, cache)
        mmin = mg.min()
        if verbose: print(f"  it{it}: L={length(template,x):.9f} mm={mmin:+.1e}", flush=True)
        if mmin > -2e-9: break
        bad = np.where(mg < min(-1e-10, mmin*0.3))[0]
        sel=[]
        if len(bad):
            for g in np.split(bad, np.where(np.diff(bad)>5)[0]+1):
                sel.append(g[np.argmin(mg[g])])
        U = np.vstack([U, Uv[sel]])
    return x, length(template,x), float(mmin)
