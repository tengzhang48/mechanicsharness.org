"""Reference zone length, load and energy ratio against speed for lambda = 0.1, d = 0.2.

Same equations as assets/case_studies_review_checks_2026-10-23/viscoelastic_check.py
(Hui and Long Eqs. 17 and 51): a = (2/pi) acos(exp(-eps/2)), R = a^2/d, eps from the
critical-separation integral. Writes eps_curve.json next to this file.
"""
import json
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

lam, d = 0.1, 0.2


def rhs(eps, V):
    def f(t):
        xi = eps * t * t
        al2 = -np.expm1(-(eps - xi)) * np.exp(-xi) / (-np.expm1(-eps))
        al = np.sqrt(max(al2, 0.0))
        return (1 - (1 - lam) * np.exp(-xi / V)) * (np.log1p(al) - np.log1p(-al)) * 2 * eps * t
    return 2 / np.pi**2 * integrate.quad(f, 0, 1, limit=400, epsabs=1e-14, epsrel=1e-13)[0]


rows = []
for V in list(np.logspace(-4, 1, 26)) + [11, 12, 12.5, 13, 13.2, 13.3, 13.38]:
    eps = optimize.brentq(lambda e: rhs(e, V) - d, 1e-6, 80, xtol=1e-13)
    a = 2 / np.pi * np.arccos(np.exp(-eps / 2))
    rows.append({"V": float(V), "epsilon": eps, "a": a, "R": a * a / d})
Path(__file__).with_name("eps_curve.json").write_text(json.dumps({"lambda": lam, "d": d, "rows": rows}, indent=1) + "\n")
print(rows[0], rows[-1])
