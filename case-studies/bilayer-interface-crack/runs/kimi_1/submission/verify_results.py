"""Reproduce the numerical checks cited in method.md."""
import cmath
import json
import math

import numpy as np
from scipy.integrate import quad

import bilayer_fe as fe


def fit(r, kline, rmax, rmin=0.001, family=2):
    mask = (r <= rmax) & (r > rmin)
    rr = r[mask, None]
    bases = [np.ones_like(rr), rr ** (1 - fe.LAMBDA), rr]
    if family >= 3:
        bases += [rr ** (2 - fe.LAMBDA), rr ** 2]
    a = np.concatenate(bases, axis=1)
    return np.linalg.lstsq(a, kline[mask], rcond=None)[0][0]


def describe(k):
    return f"|K|={abs(k):.9f}, psi={math.degrees(cmath.phase(k)):.6f} deg"


with open("results.json") as f:
    ans = json.load(f)

# Independent numerical quadrature of the far-field beam energies.
eb1 = fe.E1 / (1 - fe.NU**2)
eb2 = fe.E2 / (1 - fe.NU**2)
ga_a = (quad(lambda y: 1 / (2*eb1), 0, 1)[0] +
        quad(lambda y: (1/fe.H + 12*((fe.H+1)/2)*(y+fe.H/2)/fe.H**3)**2 /
             (2*eb2), -fe.H, 0)[0])
ga_b = (quad(lambda y: (12*(y-.5))**2/(2*eb1), 0, 1)[0] +
        quad(lambda y: (12*(y+fe.H/2)/fe.H**3)**2/(2*eb2), -fe.H, 0)[0])
assert math.isclose(ga_a, ans['caseA']['G_E2_h_over_P2'], rel_tol=1e-13)
assert math.isclose(ga_b, ans['caseB']['G_E2_h3_over_M2'], rel_tol=1e-13)
assert math.isclose(fe.G_FACTOR, 64/135, rel_tol=1e-13)
for name, g in [('caseA', ga_a), ('caseB', ga_b)]:
    kmag = math.sqrt(g / fe.G_FACTOR)
    key = ('Kabs_sqrt_h_over_P' if name == 'caseA' else 'Kabs_h32_over_M')
    assert math.isclose(kmag, ans[name][key], rel_tol=1e-13)

# Superposition check for the combined loading M=Ph=P (h=1).
ka = ans['caseA']['Kabs_sqrt_h_over_P'] * cmath.exp(
    1j*math.radians(ans['caseA']['psi_deg']))
kb = ans['caseB']['Kabs_h32_over_M'] * cmath.exp(
    1j*math.radians(ans['caseB']['psi_deg']))
assert math.isclose(math.degrees(cmath.phase(ka+kb)), ans['psi_both_deg'],
                    abs_tol=1e-12)

print('Analytic energy, K/G identity, and superposition checks: PASS')

fine = np.load('run_fine.npz')
for family, rmax in [(2, .16), (2, .32), (2, .64), (3, .16), (3, .32)]:
    k = fit(fine['r'], fine['Kline'], rmax, family=family)
    print(f'fine openings, family={family}, rmax={rmax}:',
          describe(k[0]), ';', describe(k[1]))
    for j, case in enumerate(('caseA', 'caseB')):
        key = ('Kabs_sqrt_h_over_P' if case == 'caseA' else 'Kabs_h32_over_M')
        assert abs(abs(k[j])/ans[case][key]-1) < .006

traction = np.load('traction_medium.npz')
for rmax in (.08, .16):
    k = fit(traction['r'], traction['K'], rmax, family=2)
    print(f'medium interface traction, rmax={rmax}:',
          describe(k[0]), ';', describe(k[1]))
    for j, case in enumerate(('caseA', 'caseB')):
        key = ('Kabs_sqrt_h_over_P' if case == 'caseA' else 'Kabs_h32_over_M')
        assert abs(abs(k[j])/ans[case][key]-1) < .006

# Truncation sensitivity: same fits for the L=6 and L=12 medium runs.
medium = np.load('run_medium.npz')
short = np.load('run_L6.npz')
k12 = fit(medium['r'], medium['Kline'], .32)
k6 = fit(short['r'], short['Kline'], .32)
assert np.max(np.abs(k12-k6)) < 1e-5
print('L=6 versus L=12 maximum complex-K fit difference:',
      np.max(np.abs(k12-k6)))
print('All saved-data verification checks: PASS')
