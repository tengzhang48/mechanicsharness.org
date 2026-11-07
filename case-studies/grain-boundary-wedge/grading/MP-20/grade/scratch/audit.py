"""Read-only numerical audit of the delivered solver on the scratch copy."""
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent / "submission"))
import solve
import gb
import gb2
import gb5
import gb7

xq, wq = leggauss(2000)
xq, wq = (xq + 1) / 2, wq / 2

root = Path(__file__).resolve().parent.parent
original = json.loads((root / "submission/results.json").read_text())
rerun = json.loads((root / "scratch/submission/results.json").read_text())
print("Exact original/rerun results.json match:", original == rerun)

print("\nSpectrum, interface curvature, and polynomial derivative BC residuals")
for n in (10, 18, 26, 30):
    nodes = np.arange(n + 1) * (1.0 / n)
    E = solve.elasticity_matrix(nodes, xq, wq)
    D2 = solve.d2_uniform(n)
    rate = -(E @ D2)[1:, 1:]
    ev, V = np.linalg.eig(rate)
    print("n", n, "min real", ev.real.min(), "max |imag|", abs(ev.imag).max(),
          "negative-real count", sum(ev.real < 0))
    pairs = sorted((ev[j].real, j) for j in range(n)
                   if abs(ev[j].imag) < 1e-3 * abs(ev[j]) and ev[j].real > 1e-3)
    v = np.r_[0.0, (V[:, pairs[0][1]] / V[-1, pairs[0][1]]).real]
    D1poly, D2poly = solve.nodal_derivatives(nodes)
    print(" lambda1", pairs[0][0], "D2_fd f(1)", (D2 @ v)[-1],
          "D1_polynomial f(1)", (D1poly @ v)[-1],
          "D2_polynomial f(1)", (D2poly @ v)[-1])
    if n == 26:
        l, j = min((ev[j].real, j) for j in range(n))
        v_bad = V[:, j].real
        print(" purported gauge eigenvalue", l,
              "initial-state coefficient", np.linalg.solve(V, np.ones(n))[j],
              "constant-direction overlap", abs(v_bad @ np.ones(n)) /
              (np.linalg.norm(v_bad) * np.sqrt(n)))
        print(" polynomial basis max derivative magnitude at quadrature",
              abs(solve.lagrange_basis(nodes)(xq)).max())

print("\nRoundoff sensitivity of algebraically equivalent node definitions")
for n in (26, 30):
    for definition in ("multiply", "divide"):
        nodes = (np.arange(n+1) * (1.0/n) if definition == "multiply"
                 else np.arange(n+1)/n)
        E = solve.elasticity_matrix(nodes, xq, wq)
        ev = np.linalg.eigvals(-(E @ solve.d2_uniform(n))[1:, 1:])
        pos = sorted(v.real for v in ev if abs(v.imag) < 1e-3*abs(v) and v.real > 1e-3)
        print("n", n, definition, "lambda1", pos[0])

print("\nEquilibrium tip fit sensitivity")
for n in (16, 24, 30, 40, 60):
    sg, _ = leggauss(n)
    nodes = (sg + 1) / 2
    E = solve.elasticity_matrix(nodes, xq, wq)
    w, *_ = np.linalg.lstsq(np.vstack([E, np.ones(n)]),
                           np.r_[np.ones(n), 0.0], rcond=None)
    print("n", n, "max elasticity residual", abs(E @ w - 1).max())
    for cutoff in (0.85, 0.95, 0.98):
        sel = nodes > cutoff
        if sum(sel) < 2:
            continue
        X = np.column_stack([np.ones(sum(sel)), np.sqrt(1 - nodes[sel])])
        coef, *_ = np.linalg.lstsq(X, w[sel], rcond=None)
        print(" cutoff", cutoff, "intercept", coef[0], "c", coef[1],
              "claimed K", coef[1] * np.pi * np.sqrt(2),
              "fit max residual", abs(X @ coef - w[sel]).max())

print("\nSymbolic auxiliary checks")
a, b = sp.symbols("a b", positive=True)
S = 1 / (a-b) - 1 / (a+b) - 2*b*(a-b)/(a+b)**3
S2 = sp.diff(S, a, 2)
print("S_aa:", sp.simplify(S2))
print("regular diagonal S_aa:", sp.simplify(S2-2/(a-b)**3).subs(a, b).simplify())
print("S(a,b)-S(b,a) at (.25,.75):",
      float(S.subs({a:.25,b:.75}) - S.subs({a:.75,b:.25})))
print("S_bb - S_aa:", sp.factor(sp.diff(S,b,2)-S2))
print("J2 code at .25:", gb2.J2(.25),
      "correct finite part:", -1/.25-1/(1-.25),
      "correct linear-subtraction coefficient:", 1/.25+1/(1-.25))
first_derive = 2/(a-b)**3-2/(a+b)**3+12*b*(a-b)/(a+b)**5-12*b/(a+b)**4
derive_target = 2/(a-b)**3-2/(a+b)**3-12*b/(a+b)**4+12*b/(a+b)**5-24*b**2/(a+b)**6
print("derive.py wrong derivative discrepancy:", sp.factor(first_derive-S2))
print("derive.py expansion discrepancy:", sp.factor(first_derive-derive_target))

print("\nWeak/Galerkin claimed symmetry")
nodes, weights = leggauss(8)
nodes = (nodes + 1)/2
W = gb5.weight_matrix(nodes, xq, wq)
print("gb5 W relative asymmetry:", np.linalg.norm(W-W.T)/np.linalg.norm(W))
K5 = gb5.weak_operator(nodes, xq, wq)
print("gb5 weak operator relative asymmetry:", np.linalg.norm(K5-K5.T)/np.linalg.norm(K5))
print("interior basis endpoint values are not zero; constant interpolation at 0 = 1")
print("gb7 returns only K, not (M,K):", type(gb7.galerkin_operator(nodes)),
      gb7.galerkin_operator(nodes).shape)
