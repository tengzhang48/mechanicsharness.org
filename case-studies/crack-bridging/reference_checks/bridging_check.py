"""Checks for the crack-bridging case page (probes T1 and T2 of P038, Bao and Suo 1992), 25 October 2026.

Recomputes the reference numbers from the stated formulas, by routes written for this page:
  T1: the finite-crack relations of Bao and Suo (1992), Eqs. 45-46 (p. 362), with K_tip = sqrt(E' Gamma0)
      and the root opening equal to delta_c, solved by bisection;
  T2: the exact beam profile and its characteristic equation, with the moment from the energy balance.
Also evaluates the comparison numbers quoted on the page, reads the run times, scores and verification points from
the bundled copies, and counts reasoning events in the four DeepSeek session logs (read at their original path;
the logs are not bundled). Standard library only. Run from this folder:  python3 -I bridging_check.py
These checks evaluate stated formulas and saved data; they are not new derivations.
"""
import json
import math
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

BUNDLE = Path(__file__).resolve().parents[2]
PROBES = BUNDLE / "evidence/other/corpus-runs-cpu-vm/probes"
LOGS = BUNDLE / "evidence/other/corpus-runs-cpu-vm/probes-logs"
SCORES = BUNDLE / "<study-folder>/corpus/probes/scores_P038.json"
ORIG_LOGS = Path("<probe-runs>-logs")

print(f"# bridging_check.py run {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ} with python {platform.python_version()}")
print("# command: python3 -I bridging_check.py (from assets/solution-evidence-2026-10-25)")


def bisect(f, lo, hi, tol=1e-15):
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < tol * max(1.0, abs(mid)):
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- 1. T1, finite bridged slit
E, nu, a, sb, dc, G0 = 130e9, 0.28, 2.40e-3, 160e6, 9.00e-6, 22.0
Ep = E / (1 - nu ** 2)
K0 = math.sqrt(Ep * G0)
ell = Ep * dc / sb


def sigma_tip(c, K0=K0):
    """Eq. 46 with K_tip = K0, solved for the remote stress."""
    return K0 / math.sqrt(math.pi * c) + 2 * sb / math.pi * math.acos(a / c)


def opening_eq45(c, S):
    eta = c / a
    return sb * a / Ep * (4 * math.sqrt(eta ** 2 - 1) * (S / sb - 2 / math.pi * math.acos(1 / eta))
                          + 8 / math.pi * math.log(eta))


def solve_t1(K0=K0):
    c = bisect(lambda c: opening_eq45(c, sigma_tip(c, K0)) - dc, a * (1 + 1e-12), a * 100)
    return sigma_tip(c, K0), c


S1, c1 = solve_t1()
L1 = c1 - a
print("\n## 1. T1: finite slit with constant bridging traction (Eqs. 45-46, K_tip = K0, delta(a) = delta_c)")
print(f"E' = {Ep / 1e9:.7f} GPa, K0 = sqrt(E' Gamma0) = {K0 / 1e6:.7f} MPa m^1/2")
print(f"material length E' delta_c / sigma_b = {ell * 1e3:.6f} mm; a / (E' delta_c / sigma_b) = {a / ell:.6f}; "
      f"Gamma0 / (sigma_b delta_c) = {G0 / (sb * dc):.7f}")
print(f"Sigma = {S1 / 1e6:.7f} MPa (key 133.7302861), L = {L1 * 1e3:.7f} mm per end (key 4.0856091), "
      f"c = {c1 * 1e3:.6f} mm, L/a = {L1 / a:.6f}")
print(f"  relative to the key: Sigma {S1 / 1e6 / 133.730286053811 - 1:.1e}, L {L1 * 1e3 / 4.08560913824553 - 1:.1e}")
print(f"initiation (L -> 0): Sigma = K0/sqrt(pi a) = {K0 / math.sqrt(math.pi * a) / 1e6:.4f} MPa")
grid = [a * (1 + 1e-6) + i * (L1 - 1e-6 * a) / 4000 for i in range(4001)]
mono_S = all(sigma_tip(grid[i + 1]) > sigma_tip(grid[i]) for i in range(4000))
mono_d = all(opening_eq45(grid[i + 1], sigma_tip(grid[i + 1])) > opening_eq45(grid[i], sigma_tip(grid[i])) for i in range(4000))
print(f"Sigma increasing up to the event: {mono_S}; root opening increasing up to the event: {mono_d}")
S0, c0 = solve_t1(K0=0.0)
u0 = math.pi * Ep * dc / (8 * sb * a)
print(f"Gamma0 = 0 (the paper's Eq. 47): Sigma = {S0 / 1e6:.4f} MPa, L = {(c0 - a) * 1e3:.4f} mm; "
      f"closed form (2/pi) acos(exp(-pi E' dc/(8 sb a))) sigma_b = {2 / math.pi * math.acos(math.exp(-u0)) * sb / 1e6:.4f} MPa")
print(f"  dropping Gamma0 changes Sigma by {S0 / S1 - 1:+.2%} and L by {(c0 - a) / L1 - 1:+.1%}")
S_lefm = math.sqrt(Ep * (G0 + sb * dc) / (math.pi * a))
print(f"toughness replacement, Gamma0 + sigma_b delta_c = {G0 + sb * dc:.0f} J/m^2 on the original slit: "
      f"{S_lefm / 1e6:.4f} MPa, {S_lefm / S1 - 1:+.2%} relative to the key")
g = G0 / (sb * dc)
L0_ssb = math.pi / 8 * ell * (math.sqrt(1 + g) + math.sqrt(g)) ** -2
print(f"small-scale estimate of the bridged length, Eq. 20: L0 = (pi/8) ell [sqrt(1+g)+sqrt(g)]^-2 = {L0_ssb * 1e3:.4f} mm "
      f"= {L0_ssb / a:.3f} a (Gamma0 = 0: (pi/8) ell = {math.pi / 8 * ell * 1e3:.4f} mm)")
# DeepSeek MP without the unit (R03): its root, and the opening Eq. 45 gives there
S_ds, L_ds = 158.524900100573e6, 34.324303918593e-3
c_ds = a + L_ds
print(f"P038-T1-without-dsaz root Sigma = {S_ds / 1e6:.6f} MPa, L = {L_ds * 1e3:.6f} mm: "
      f"its formula (4/E')[Sigma sqrt(c^2-a^2) - t (c-a)] = {4 / Ep * (S_ds * math.sqrt(c_ds ** 2 - a ** 2) - sb * L_ds) * 1e6:.6f} um; "
      f"Eq. 45 = {opening_eq45(c_ds, S_ds) * 1e6:.6f} um")
print(f"  its K_tip at that root, Eq. 46: {math.sqrt(math.pi * c_ds) * (S_ds - 2 * sb / math.pi * math.acos(a / c_ds)) / 1e6:.6f} MPa m^1/2 "
      f"(K0 = {K0 / 1e6:.6f})")
print(f"  relative errors: Sigma {S_ds / S1 - 1:+.2%}, L {L_ds / L1 - 1:+.1%}")

# ---------------------------------------------------------------- 2. T2, unequal strips
b, E1, h1, E2, h2 = 25.0e-3, 110e9, 0.60e-3, 190e9, 0.90e-3
sb2, dc2, G02 = 0.20e6, 0.20e-3, 1.40
gam = G02 / (sb2 * dc2)
kap = math.sqrt(2 * gam)


def char(lam):
    return 1 + math.cosh(lam) * math.cos(lam) - kap * math.sinh(lam) * math.sin(lam)


lam = bisect(char, 1.0, 1.8751)
lam0 = bisect(lambda z: 1 + math.cosh(z) * math.cos(z), 1.0, 2.0)
A = ((math.sinh(lam) + math.sin(lam)) - kap * (math.sinh(lam) - math.sin(lam))) / (math.cosh(lam) + math.cos(lam))


def u(z):
    return (1 - (math.cosh(z) + math.cos(z)) / 2 + kap / 2 * (math.cosh(z) - math.cos(z))
            + A / 2 * (math.sinh(z) - math.sin(z)))


def u1(z):
    return (-(math.sinh(z) - math.sin(z)) / 2 + kap / 2 * (math.sinh(z) + math.sin(z))
            + A / 2 * (math.cosh(z) - math.cos(z)))


def u3(z):
    return (-(math.sinh(z) + math.sin(z)) / 2 + kap / 2 * (math.sinh(z) - math.sin(z))
            + A / 2 * (math.cosh(z) + math.cos(z)))


def arms(D1, D2):
    C = 1 / D1 + 1 / D2
    ellb = (dc2 / (C * sb2)) ** 0.25
    m = math.sqrt(2 * (G02 + sb2 * dc2 / 2) / C)
    return {"C": C, "ell": ellb, "M": b * m, "L": lam * ellb, "dmid": dc2 * u(lam / 2), "d5": dc2 * u(5e-3 / ellb),
            "rot": dc2 / ellb * u1(lam)}


D1, D2 = E1 * h1 ** 3 / 12, E2 * h2 ** 3 / 12
key = arms(D1, D2)
Dav = ((E1 + E2) / 2) * ((h1 + h2) / 2) ** 3 / 12
avg = arms(Dav, Dav)
one = arms(D1, 1e30)
print("\n## 2. T2: unequal strips, linear softening connectors, front work")
print(f"D1 = {D1:.6f} N m, D2 = {D2:.6f} N m (per width); C = 1/D1 + 1/D2 = {key['C']:.7f} (N m)^-1")
print(f"gamma = Gamma0/(sigma_b delta_c) = {gam:.4f}, kappa = sqrt(2 gamma) = {kap:.6f}, ell = (delta_c/(C sigma_b))^1/4 = {key['ell'] * 1e3:.6f} mm")
print(f"lambda (smallest root of 1 + cosh cos - kappa sinh sin = 0) = {lam:.12f}, A = {A:.12f}; u(lambda) - 1 = {u(lam) - 1:.1e}, "
      f"u'''(lambda) = {u3(lam):.1e}")
print(f"energy: C m^2/2 = Gamma0 + sigma_b delta_c/2 = {G02 + sb2 * dc2 / 2:.2f} J/m^2 -> m = {key['M'] / b:.7f} N, M = {key['M']:.10f} N m (key 0.2126257874)")
print(f"L = {key['L'] * 1e3:.8f} mm (key 10.80010237), delta(L/2) = {key['dmid'] * 1e6:.8f} um (key 38.21889793)")
print(f"delta(x = 5 mm) = {key['d5'] * 1e6:.6f} um; rear relative rotation = {key['rot']:.5f} rad; L / h2 = {key['L'] / h2:.4f}")
prof = [u(lam * i / 400) for i in range(401)]
print(f"profile on [0, L]: min {min(prof):.3e}, max {max(prof):.6f}, increasing: {all(prof[i + 1] > prof[i] for i in range(400))}")
print(f"Gamma0 = 0: first cantilever root 1 + cosh cos = 0 at {lam0:.12f}")
print(f"averaged equal arms: M = {avg['M']:.5f} N m, L = {avg['L'] * 1e3:.4f} mm, delta(L/2) = {avg['dmid'] * 1e6:.6f} um, "
      f"delta(5 mm) = {avg['d5'] * 1e6:.4f} um")
print(f"upper-arm compliance only: M = {one['M']:.5f} N m, L = {one['L'] * 1e3:.4f} mm, delta(L/2) = {one['dmid'] * 1e6:.6f} um, "
      f"delta(5 mm) = {one['d5'] * 1e6:.4f} um")
print("  delta(L/2) = delta_c u(lambda/2) depends on gamma only; L and delta at a fixed x scale with ell ~ C^(-1/4)")
Eh = 150e9
h = 1e-3
Ceq = 2 / (Eh * h ** 3 / 12)
print(f"equal arms, D = E h^3/12 with E = 150 GPa, h = 1 mm, m = 1 N: C m^2/2 = {Ceq / 2:.6e} J/m^2; 12 m^2/(E h^3) = {12 / (Eh * h ** 3):.6e} J/m^2 "
      "(the form of the paper's Eq. 51 left side, which has the plane-strain E' where T2 has E)")

# ---------------------------------------------------------------- 3. runs: times, scores, verification points
print("\n## 3. Runs: minutes from start/exit stamps, grader totals, verification rows")
scores = json.loads(SCORES.read_text())
for name in sorted(k for k in scores if k.split("-")[1] in ("T1", "T2")):
    st = datetime.fromisoformat((LOGS / f"{name}.start").read_text().strip().replace("Z", "+00:00"))
    ex_txt = (LOGS / f"{name}.exit").read_text().split()
    ex = datetime.fromisoformat(ex_txt[2].replace("Z", "+00:00"))
    lab = scores[name]["label"]
    grade = (PROBES / "grading" / lab / "GRADE.md").read_text()
    ver = re.search(r"\|\s*Verification\s*\|\s*([\d.]+)\s*/\s*(\d+)", grade)
    print(f"{name:24s} {lab}  exit {ex_txt[1]}  {(ex - st).total_seconds() / 60:5.1f} min  total {scores[name]['total']:5.1f}  "
          f"verification {ver.group(1)}/{ver.group(2)}")

# ---------------------------------------------------------------- 4. reasoning events in the DeepSeek logs
print("\n## 4. Reasoning events in the DeepSeek session logs (original path; logs not bundled)")
for name in ("P038-T1-with-dsaz", "P038-T1-without-dsaz", "P038-T2-with-dsaz", "P038-T2-without-dsaz"):
    p = ORIG_LOGS / f"{name}.jsonl"
    if not p.exists():
        print(f"{name}: log not found at {p}")
        continue
    txt = p.read_text()
    print(f"{name}: {txt.count(chr(34) + 'contentType' + chr(34) + ':' + chr(34) + 'reasoning' + chr(34))} reasoning events, "
          f"{txt.count(chr(34) + 'contentType' + chr(34) + ':' + chr(34) + 'text' + chr(34))} text events")
