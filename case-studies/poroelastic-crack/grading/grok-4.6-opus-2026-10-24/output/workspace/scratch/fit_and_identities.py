"""Check polynomial fits and a few analytic identities. Read-only on C/."""
import json
import os
import numpy as np
from scipy.special import gamma as Gf

HERE = os.path.dirname(os.path.abspath(__file__))

rows = json.load(open(os.path.join(HERE, "a_chat_vs_gamma.json")))
fx = json.load(open(os.path.join(HERE, "final_extrapolation.json")))
rows.append({"gamma": 0.25, "a": fx["a"], "chat": fx["chat"]})
rows.sort(key=lambda r: r["gamma"])

print("gamma    a_num     a_fit     da        chat_num  chat_fit  dchat")
max_da = 0.0
max_dc = 0.0
for r in rows:
    g = r["gamma"]
    af = 1 - 0.75 * g - 0.0634 * g**2 - 0.0582 * g**3
    cf = 2 * np.pi * g + 0.434 * g**2 + 0.37 * g**3
    da = af - r["a"]
    dc = cf - r["chat"]
    max_da = max(max_da, abs(da))
    max_dc = max(max_dc, abs(dc))
    print(f"{g:6.3f}  {r['a']:.7f}  {af:.7f}  {da:+.2e}  {r['chat']:.5f}  {cf:.5f}  {dc:+.2e}")
print(f"max |da|={max_da:.3e}  claimed <1e-4")
print(f"max |dchat|={max_dc:.3e}  claimed <5e-4")

print("\nbounds at tabulated gamma:")
for r in rows:
    g = r["gamma"]
    lo = (1 - g) ** 0.75
    hi = np.sqrt((1 - g) * (1 - g / 2))
    a = r["a"]
    print(f"g={g:.3f}  lo={lo:.6f}  a={a:.6f}  hi={hi:.6f}  in=[{lo < a < hi}]")

print("\nYu-like (1-gamma)^0.743 vs a")
for r in rows:
    g = r["gamma"]
    yu = (1 - g) ** 0.743
    print(f"g={g:.3f}  yu={yu:.6f}  a={r['a']:.6f}  rel={(yu - r['a']) / r['a']:+.4f}")

nu, nuu, al = 0.2, 0.4, 0.89
eta = al * (1 - 2 * nu) / (2 * (1 - nu))
print("\n4 ln 2 - 2 gamma_E =", 4 * np.log(2) - 2 * np.euler_gamma)
print("V2 small-time coeff =", -16 / np.sqrt(np.pi) * eta * (1 - nuu))
print("Laplace lim sqrt(s)*s V2 =", -8 * eta * (1 - nuu))
print("K2_t14 from json", fx["K2_t14_coeff"])
print("Gamma(5/4)=", Gf(1.25), " claimed s^{1/4}sK2 -> -0.23554")
print("time coeff from -0.23554 / Gamma(1.25) =", -0.23554 / Gf(1.25))
