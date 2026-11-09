"""Reproduce TOPIC §3 existence table from williams_exponent."""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent / "crack_tip_asymptotics" / "checks"))
from crack_tip_checks import williams_exponent

crit2b = brentq(lambda x: np.tan(x) - x, 4.0, 4.6)
print(f"tan(2beta)=2beta at 2beta={np.degrees(crit2b):.5f} deg")

rows = [180, 240, 257.45, 270, 300, 360]
print(f"{'2beta':>10} {'pI':>12} {'pII':>14}")
for ang in rows:
    beta = np.pi if ang == 360 else np.radians(ang) / 2.0
    pI = williams_exponent(beta, "I")
    try:
        pII = williams_exponent(beta, "II")
        pIIs = f"{pII:.6f}"
    except Exception as e:
        pIIs = f"CRASH {e}"
    print(f"{ang:10} {pI:12.6f} {pIIs:>14}")

print("exact crit mode II:", end=" ")
try:
    print(williams_exponent(crit2b / 2, "II"))
except Exception as e:
    print(type(e).__name__, e)
