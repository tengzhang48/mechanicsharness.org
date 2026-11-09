"""check_field at the table's mode-II threshold angle."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "crack_tip_asymptotics" / "checks"))
from crack_tip_checks import check_field, williams_exponent

beta = np.radians(257.45) / 2
th = np.linspace(0, beta, 801)
# dummy even/odd-ish field, p=0
srr = np.zeros_like(th)
stt = np.zeros_like(th)
srt = np.cos(np.pi * th / (2 * beta))  # srt even-ish for mode II
print("williams pII", williams_exponent(beta, "II"))
try:
    rep = check_field(th, srr, stt, srt, p=0.0, law="linear", beta=beta, mode="II")
    print("check_field returned", len(rep), "items")
except Exception as e:
    print("check_field CRASH", type(e).__name__, e)

# 90 mode II p=-1 field? skip
print("90 II p", williams_exponent(np.pi/2, "II"))
# A rigid rotation is excluded; p=-1 is lambda=2.
# Construct mode II Airy at beta=pi/2, lam=2, p=-1
beta = np.pi / 2
p = williams_exponent(beta, "II")
lam = 1 - p
print("lam", lam, "p", p)
s = lam + 1
th = np.linspace(0, beta, 2001)
B = -np.sin((lam + 1) * beta) / np.sin((lam - 1) * beta)
print("B", B)
F = np.sin((lam + 1) * th) + B * np.sin((lam - 1) * th)
Fp = (lam + 1) * np.cos((lam + 1) * th) + B * (lam - 1) * np.cos((lam - 1) * th)
Fpp = -(lam + 1)**2 * np.sin((lam + 1) * th) - B * (lam - 1)**2 * np.sin((lam - 1) * th)
try:
    from crack_tip_checks import format_report
    rep = check_field(th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p=p, law="linear",
                      beta=beta, mode="II")
    fails = {r["check"] for r in rep if r["passed"] is False}
    print("90 II field failed", fails)
    for r in rep:
        if "exponent" in r["check"] or "non-singular" in r["check"] or "leading" in r["check"]:
            print(r["check"], r["passed"], r.get("note"), r["value"])
except Exception as e:
    print("90 II field CRASH", type(e).__name__, e)
