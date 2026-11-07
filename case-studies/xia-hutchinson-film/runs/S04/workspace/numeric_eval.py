"""Numerical evaluation of the requested normalized quantities.

Normalized quantities (analytic results, verified symbolically in analytic_checks.py):
  gamma_A = G_A Ebar/(sigma0^2 ell)   = tanh(H/(2ell))
  gamma_B = G_B Ebar/(sigma0^2 ell)   = 2 tanh(H/(2ell)) - tanh(H/ell) = 2 t^3/(1+t^2)
  kI_A    = KI_A/(sigma0 sqrt(ell))   = sqrt((1-nu^2) tanh(H/(2ell)))       (KII = 0 by symmetry)
  kI_B    = KI_B/(sigma0 sqrt(ell))   = sqrt((1-nu^2) [2 tanh(H/(2ell)) - tanh(H/ell)])
Sampling at ell/H = 0.5, 1, 1.5, 2  =>  H/(2ell) = 1/(2*(ell/H)) = 1, 0.5, 1/3, 0.25.
nu = 0.3 as requested.
"""
import json
import numpy as np

NU = 0.3
rows = []
for eoH in (0.5, 1.0, 1.5, 2.0):
    x = 1.0 / (2.0 * eoH)          # H/(2 ell)
    t = np.tanh(x)
    gA = t
    gB = 2.0 * t - np.tanh(2.0 * x)
    gB_alt = 2.0 * t**3 / (1.0 + t**2)
    assert abs(gB - gB_alt) < 1e-15
    kA = np.sqrt((1 - NU**2) * gA)
    kB = np.sqrt((1 - NU**2) * gB)
    rows.append(dict(ell_over_H=eoH, H_over_2ell=x,
                     gamma_A=gA, gamma_B=gB, kI_A=kA, kI_B=kB))

# critical stress: sigma0c = sqrt(Ebar Gc / (ell * gamma)) ; report the dimensionless
# combination sigma0c sqrt(ell/(Ebar Gc)) = 1/sqrt(gamma)
print('ell/H   H/2ell    gamma_A    gamma_B    kI_A       kI_B       sigcA*     sigcB*')
for r in rows:
    scA = 1.0 / np.sqrt(r['gamma_A'])
    scB = 1.0 / np.sqrt(r['gamma_B'])
    print(f"{r['ell_over_H']:4.1f}  {r['H_over_2ell']:7.4f}  {r['gamma_A']:9.6f}  {r['gamma_B']:9.6f}  "
          f"{r['kI_A']:9.6f}  {r['kI_B']:9.6f}  {scA:9.6f}  {scB:9.6f}")
print('* sigc = sigma0c sqrt(ell/(Ebar Gc)) = gamma^{-1/2}')

with open('numeric_rows.json', 'w') as f:
    json.dump(rows, f, indent=2)
print('wrote numeric_rows.json')
