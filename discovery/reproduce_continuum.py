"""Regenerates the continuum smoothed-staple optimum at beta = 36 deg
(the L_infinity = 1.2826726243 value quoted in the paper), addressing the
reproducibility gap that results/arc2_final.npy had no pipeline.

Runtime: ~10-20 min. This is DISCOVERY code (floating point), not a proof.
Expected final output: L ~= 1.2826726 with dense support margin ~ -1e-14.
"""
import numpy as np
from arcmodel2 import solve, build, tri

d = np.deg2rad
beta = d(36.0)
T = ['s','a','s','a','c','s','c','a','s','a','s']
seed = [d(-110), 0.05, 0.95, d(6), 0.05, 0.59, d(6), d(77), 0.36,
        d(77), 0.59, d(6), 0.05, 0.95, d(6), 0.05]
x, L, mm = solve(T, seed, beta, rounds=14, m0=720, verbose=True)
# continue with a denser grid for the final digits
x, L, mm = solve(T, x, beta, rounds=10, m0=5760, verbose=True)
print(f"\nfinal: L = {L:.10f}   dense support minmargin {mm:+.2e}")
names = ["theta0","s1","r1a","phi1a","s1m","r1b","phi1b","delta1","base",
         "delta2","r2b","phi2b","s2m","r2a","phi2a","s3"]
for nm, v in zip(names, x):
    if nm == "theta0" or "phi" in nm or "delta" in nm:
        print(f"  {nm:7s} = {np.degrees(v):+9.4f} deg")
    else:
        print(f"  {nm:7s} = {v:9.6f}")
print("\nreference: sin(36 deg) =", np.sin(beta))
np.save('results/arc2_final_reproduced.npy', x)
