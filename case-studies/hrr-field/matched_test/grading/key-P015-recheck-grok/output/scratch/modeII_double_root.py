"""How close to the mode-II double root does williams_exponent crash?"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent / "crack_tip_asymptotics" / "checks"))
from crack_tip_checks import williams_exponent, check_field

crit2b = brentq(lambda x: np.tan(x) - x, 4.0, 4.6)
print("crit 2beta deg", np.degrees(crit2b))

for ddeg in [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.0]:
    for sign in (-1, 1):
        ang = np.degrees(crit2b) + sign * ddeg
        beta = np.radians(ang) / 2
        try:
            p = williams_exponent(beta, "II")
            print(f"2beta={ang:.6f} p={p:.8f}")
        except Exception as e:
            print(f"2beta={ang:.6f} CRASH {e}")

# test_checks offsets: crit is beta, +-0.005 rad
crit_beta = crit2b / 2
for db in (-0.005, 0.005, -0.001, 0.001, -0.0001, 0.0001):
    try:
        p = williams_exponent(crit_beta + db, "II")
        print(f"beta=crit{db:+g} 2beta_deg={np.degrees(2*(crit_beta+db)):.5f} p={p:.8f}")
    except Exception as e:
        print(f"beta=crit{db:+g} CRASH {e}")

# T-stress opening robustness: many random grids?
# skip

# Does leading T-stress FAIL only crack exponent? yes.

# Opening: T-stress with explicit 0 displacement should fail? hydrostatic Ut=0 failed.
# For a crack T-stress, is opening physically required? No. The check still applies
# because crack=True and mode I. Roundoff made it pass.

# Check g(lam) near double root to see exclusion of lam=1
beta = crit2b / 2
g = lambda lam: np.sin(2 * lam * beta) - lam * np.sin(2 * beta)
print("g(1) at crit", g(1.0))
print("g'(via fd)", (g(1+1e-8)-g(1-1e-8))/2e-8)

# at 270 mode II p=0.09147 table 0.092
print("270 II", williams_exponent(np.radians(135), "II"))
print("240 II", williams_exponent(np.radians(120), "II"))  # none singular
