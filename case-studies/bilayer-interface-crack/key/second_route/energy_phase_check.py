"""Exact check of the RELATIVE phase psi_A - psi_B from beam-theory energy alone.
K is linear in (P, M): K = P K_A + M K_B (h = 1), and G(P, M) = Gfac |K|^2 = a P^2 + 2 b P M + c M^2 from the beam
expression.  Hence cos(psi_A - psi_B) = Re(K_A conj K_B)/(|K_A||K_B|) = b / sqrt(a c), independent of any FE.
Also derives the load ratio M/(P h) at which psi = 0 (if it exists for P > 0)."""
import json
import numpy as np
from bilayer import SETS

fin = json.load(open("results/final.json"))
out = {}
print("| set | psi_A - psi_B (FE, M-int) | arccos(b/sqrt(ac)) (beam energy, exact) | difference | M/(P h) at psi = 0 |")
print("|---|---|---|---|---|")
for s in ["V1", "S1", "S2", "S3", "S4"]:
    E1, nu1, E2, nu2, hH = SETS[s]
    h, H = 1.0, 1.0 / hH
    Eb1, Eb2 = E1 / (1 - nu1 ** 2), E2 / (1 - nu2 ** 2)
    a = (1 / h) / (2 * Eb1) + (1 / H + 3 * (h + H) ** 2 / H ** 3) / (2 * Eb2)
    c = 12 / (2 * Eb1 * h ** 3) + 12 / (2 * Eb2 * H ** 3)
    b = 3 * (h + H) / (Eb2 * H ** 3)
    dpsi_exact = np.degrees(np.arccos(b / np.sqrt(a * c)))
    KA = complex(*fin[s]["A"]["K"]); KB = complex(*fin[s]["B"]["K"])
    dpsi_fe = fin[s]["A"]["psi"] - fin[s]["B"]["psi"]
    chi = -KA.imag / KB.imag if abs(KB.imag) > 1e-9 else np.nan
    if np.isfinite(chi):
        Kc = KA + chi * KB
        chi_txt = f"{chi:.5f}" if (chi > 0 and Kc.real > 0) else f"none for P>0, M>0 (Im K = 0 at M/(Ph) = {chi:.4f}, where Re K = {Kc.real:.3f} < 0 or M < 0)"
    else:
        chi_txt = "0 (psi_B = 0: only M = 0 ... psi = 0 needs P = 0)"
    out[s] = dict(dpsi_fe=dpsi_fe, dpsi_exact=dpsi_exact, chi=chi)
    print(f"| {s} | {dpsi_fe:.4f} | {dpsi_exact:.4f} | {dpsi_fe - dpsi_exact:+.4f} | {chi_txt} |")
json.dump(out, open("results/energy_phase_check.json", "w"), indent=1)
