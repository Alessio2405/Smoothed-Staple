# The smoothed staple: improved certified escape paths for Bellman's lost-in-a-forest problem in isosceles triangles

Machine-checkable proofs, discovery code, and a draft manuscript for a strict
improvement of the *staple* upper bounds of

> A. Temerev, *An improved upper bound for Bellman's lost-in-a-forest problem
> in isosceles triangles* (2026), repository
> [`atemerev/bellman-staple`](https://github.com/atemerev/bellman-staple).

This work builds directly on that repository's exact reduction and
certification methodology, and improves its results.

> **Status: work in progress — not peer-reviewed.**
> What is *proved* here is exactly one thing: a certified upper bound
> ℓ(T) ≤ 1.2826879687 at the golden gnomon, improving the previous certified
> bound 1.2849615334, via explicit machine-checkable certificates.
> The smoothed line–arc geometry behind the witness is due to Gibbs (2016) —
> see *Prior work* below. Everything else in this repository (the continuum
> value, the sin β curvature law, the angle sweeps) is numerical evidence or
> conjecture, clearly labelled as such. The manuscript in `paper/` is an
> informal working draft circulated for verification only. Scrutiny and
> corrections are welcome via issues.

## What this is

For the isosceles triangle T(β) with unit legs and base angle β, an escape
path is a curve guaranteed to reach ∂T from any starting point and heading.
Escape depends only on the convex hull H of the path and reduces — via the
planar case of Lutwak's containment theorem, exactly as in `bellman-staple` —
to one polygon-contains-disk inclusion: disk(ρ) ⊆ W, where
W = 2cos β·R₀(H) + R₁(H) + R₂(H) is a weighted Minkowski sum over the edge
normals of T and ρ = 2·Area = sin 2β.

At the golden gnomon, the three-segment staple is not optimal among
convex-position polygonal paths satisfying the same escape criterion: free
polylines are strictly shorter. Numerically, the limiting shape replaces the staple's straight
shoulders by circular arcs — a line–arc geometry proposed by Gibbs in 2016
(see *Prior work*) — with observed curvature radius sin β (conjectural).

## Main results

**Theorem (certified).** At the golden gnomon T (apex 108°, base angles 36°):

```
ℓ(T) ≤ 1.2826879686692367…  <  1.2849615334…  (staple bound, Temerev 2026)
```

Witness: an explicit 20-vertex path with rational coordinates (denominator
10⁸), given in `verification/certificate20.json`. The certificate consists of
57 edge conditions (cycle convexity, positive orientation, and
D = cross² − ρ²|e|² > 0) for a convex polygon whose vertices lie in W; every
accepted sign determination is exact in **Q(√5, √(10−2√5))** via sympy. Since
the polygon's vertices are points of W and W is convex, the polygon's hull is
contained in W, so disk(ρ) ⊆ W. No floating point enters any accepted check: the length comparison uses
exact rational majorants of each square root (integer arithmetic), and
decimal printouts are display only.

A simpler, fully self-contained 8-vertex certificate
(`verification/smoothed_staple_certificate.py`, ℓ(T) ≤ 1.2827871) is included
for easy independent auditing.

**Structure (numerical, conjectured exact).** Free-radius line/arc
optimization converges, at β = 36°, to the symmetric path

```
straight 0.291027 → arc, fitted radius 0.587785 ≈ sin 36°, turning 9.5232° → straight 0.071683
  → corner 76.3054° → base 0.361665 → (mirror image),     L∞ = 1.2826726243
```

with dense support-margin −6·10⁻¹⁵. The arc radius matches sin 36° =
0.5877853 to the optimizer's precision, with the radius left free.

**Curvature law (conjecture, with an LP argument).** The escape constraint is
linear in the support function h_H, so the length-minimization is an
infinite-dimensional linear program. On active direction windows,
h_W + h_W″ = ρ forces ∂W to be a ρ-arc, i.e. Σ wᵢ·(curvature radius of ∂H at
the rotated angle) = ρ; the bang-bang optimum concentrates curvature on the
copy with the largest weight w₀ = 2cos β, predicting curvature radius
**ρ/(2cos β) = sin β** on all active arcs of the optimal escape path.

**Interval widening (numerical, conservative).** Per-angle values obtained by
grid optimization of 8-vertex polylines followed by a conservative scaling
correction, evaluated in floating point against the hull-edge criterion
(numerical, not symbolic — unlike the certificates above):

| quantity | value |
|---|---|
| gnomon certified bound (this repo) | ℓ(T) ≤ 1.2826879687 |
| previous certified staple bound | 1.2849615334 |
| continuum smoothed-staple value at 36° (numerical) | 1.2826726243 |
| gain over the per-angle optimal staple, β ∈ [30°, 33°] | ≈ +0.0026 |
| smoothed-vs-caliper crossing (conservative) | β ≈ 30.76° |
| staple-vs-caliper crossing (Temerev) | β ≈ 31.05° |
| smoothed-vs-zee crossing (resolved, n=10) | β ≈ 42.284° (staple: 42.28°) |

The smoothed family also beats the low-regime square path at every sampled
β down to 29°. On the high end the smoothing gain genuinely collapses (+1.6·10⁻⁴ at 42.0°,
+1·10⁻⁵ at 42.5°; resolved with n=10 polylines and conservative numerical correction,
`discovery/high_resolve.py`): the arcs switch off as β approaches the zee
crossing, which stays at β ≈ 42.284° — essentially the staple's 42.28°, short
of Gibbs's conjectured 42.3°. Square-path comparisons are only made
below Ward's regime switch β₀ ≈ 39.13°.

**Optimality is not claimed** for any construction here.

## Prior work and what is actually new

Philip Gibbs proposed already in 2016, on Monte Carlo evidence, non-polygonal
"Tunnel" solutions with line–arc structure for base angles roughly
30.8°–38.1° (a range containing 36°), and derived coupled curvature
conditions for extremal escape paths in convex polygons. **The qualitative
smoothed line–arc geometry is therefore his, not ours.** Those proposals were
conjectural and not accompanied by certified bounds; Temerev (2026) later
gave the first certified interval bounds via the piecewise-linear staple.
What is new here is: (a) the first exact rational polygonal certificate
showing that a smoothed/Tunnel-type path *rigorously* improves the certified
staple bound at the golden gnomon; and (b) the specific quantitative
prediction r = ρ/(2cos β) = sin β for the arc curvature radius, stated as a
conjecture.

## Repository layout

```
paper/          smoothed-staple.tex — informal working draft (amsart),
                smoothed-staple.pdf — compiled version of the same draft.
                circulated for verification only; not submitted anywhere.
                It exists to state precisely what the certificates verify
                and what they assume, so the argument can be checked.
verification/   machine-checkable certificates
  verify_certificate20.py       verifies the 20-vertex gnomon theorem (~15 min)
  certificate20.json            exact rational vertices + 40-digit length
  smoothed_staple_certificate.py  self-contained 8-vertex certificate (~2 min)
  certificate8.json             its data
discovery/      numerical code that found the witnesses (not used in proofs)
                note: sweep_high2.py and high_resolve.py import staple_family
                from a sibling clone of atemerev/bellman-staple (runtime
                dependency only; no code from that repository is included here)
  explore.py                    polyline optimizer against exact hull margins
  fast.py                       support-function grid relaxation
  arcmodel2.py                  line/arc template model + constraint generation
  certify20.py                  builds the rational certificate from P20
  sweep_fast.py, sweep_high2.py per-angle sweeps with conservative numeric correction
  reproduce_continuum.py        regenerates the continuum line/arc optimum
  results/                      solved paths, sweep data
figures/        path comparison and margin sweep figures
```

## Reproducing the proofs

Requirements: Python 3.9+, `sympy`; `numpy`/`scipy` are used only to propose
certificate cycles and in discovery code (accepted checks are symbolic
or exact integer arithmetic).

```
cd verification
python3 smoothed_staple_certificate.py   # ~2 min;  expected: ALL CERTIFICATES PASS
python3 verify_certificate20.py          # ~15 min; expected: ALL CERTIFICATES PASS
```

## Verification notes

- The certificates import exactly one external mathematical fact: the escape
  reduction disk(ρ) ⊆ W ⟹ escape, proved in Temerev (2026) via Lutwak's
  containment theorem. Everything else is verified in exact arithmetic in
  this repository.
- What is theorem-level: the two gnomon certificates. What is numerical:
  the continuum value, the sin β curvature law, and the sweep/crossing
  figures.
- Verification is deterministic: the certificate cycles (57 resp. 24 index
  triples) are stored in the JSON files, so the checkers do not depend on
  the Qhull/scipy version. A wrong stored cycle cannot pass — the exact
  convexity, orientation, and distance checks would fail on it.
- This work has not been peer-reviewed. The two places that most deserve
  scrutiny: whether the imported reduction is applied under the intended
  conventions, and the polygon-of-W-points variant of the verification
  scheme (sufficient because the hull of a subset of W lies in W).
  Issues and corrections are welcome.

## References

- A. Temerev, *An improved upper bound for Bellman's lost-in-a-forest problem
  in isosceles triangles* (2026), github.com/atemerev/bellman-staple
- P. E. Gibbs, *Lost in an Isosceles Triangle*, working paper (July 2016)
- P. E. Gibbs, *Bellman's Escape Problem for Convex Polygons*, working paper (June 2016)
- E. Lutwak, *Containment and circumscribing simplices*, Discrete Comput.
  Geom. 19 (1998), 229–235
- S. R. Finch, J. E. Wetzel, *Lost in a forest*, Amer. Math. Monthly 111
  (2004), 645–654
- J. W. Ward, *Exploring the Bellman forest problem*, manuscript (2008)

## Author and license

Alessio Doria. Code and data are released under the MIT license (see
LICENSE); the manuscript in `paper/` is © the author. No source file from bellman-staple is included verbatim. Some verification and discovery code follows or adapts its mathematical and implementation structure; attribution and the original MIT notice are provided in THIRD_PARTY_NOTICES.md. Developed with AI assistance (Claude, Anthropic).
