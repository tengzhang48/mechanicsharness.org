"""Verify J(X) = +e^{-X} Ei(X) via  PV int_{-X}^{X} (exact symmetrization) + E1 tail."""
import numpy as np
from scipy.integrate import quad
from scipy.special import expi, exp1

for X in [0.5, 1.0, 3.0]:
    sym = quad(lambda v: (np.exp(-v) - np.exp(v))/v, 0, X, limit=400)[0]
    tail = exp1(X)
    pv_total = sym + tail
    Jnum = -np.exp(-X)*pv_total
    print(f"J({X}): numeric={Jnum:.12f}  closed={np.exp(-X)*expi(X):.12f}")
