import numpy as np
from scipy.optimize import minimize
np.seterr(all='ignore')

def setup(beta, ndirs=1440):
    sb, cb = np.sin(beta), np.cos(beta)
    normals = np.array([[0,-1.],[sb,cb],[-sb,cb]])
    weights = np.array([2*cb, 1., 1.])
    rho = 2*sb*cb
    Rs = np.array([[[n[0],n[1]],[-n[1],n[0]]] for n in normals])  # R_i
    th = np.linspace(0, 2*np.pi, ndirs, endpoint=False)
    U = np.c_[np.cos(th), np.sin(th)]                              # (m,2)
    # precompute rotated directions: V_i = R_i^T u  -> (3, m, 2)
    V = np.einsum('ikl,mk->iml', Rs, U)
    return dict(V=V, w=weights, rho=rho, U=U)

def support_margins(P, S):
    # h_W(u_k) - rho for all k;  h_W(u)= sum_i w_i max_j P_j . (R_i^T u)
    # dots: (3, m, n) = V (3,m,2) . P^T (2,n)
    dots = np.einsum('imk,jk->imj', S['V'], P)
    h = (S['w'][:,None] * dots.max(axis=2)).sum(axis=0)
    return h - S['rho']

def path_length(x, n):
    P = x.reshape(n,2)
    return np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1))

def optimize(beta, n, x0, S, maxiter=800):
    cons = [dict(type='ineq', fun=lambda x: support_margins(x.reshape(n,2), S))]
    res = minimize(path_length, np.asarray(x0,float), args=(n,), constraints=cons,
                   method='SLSQP', options=dict(ftol=1e-14, maxiter=maxiter))
    P = res.x.reshape(n,2)
    return path_length(res.x, n), P, support_margins(P, S).min()
