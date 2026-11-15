"""omega(0, 0, eta) of identical layers from fe_bilayer.py: the two limits and the sensitivity to the numerical choices.

Known values: eta = 1 exact 49.107 deg (Suo and Hutchinson 1990, Eq. 2.19); eta -> 0: 52.08 deg (their computation; 52.14
from Thouless et al. 1987); Hutchinson and Suo (1992, Eq. 3.10): omega ~ 52.1 - 3 eta deg, "error believed to be within
one percent". Their 1990 Tables 2 and 3 print 54.6 (eta = 0.1) and 49.8 (eta = 0.5).
At eta = 0.5 the same case is run with a finer mesh, a longer model and two other fit windows.
Writes homogeneous.json. Run: python3 -I check_homogeneous.py
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fe_bilayer as fe  # noqa: E402
from validate_fe import gamma_deg  # noqa: E402


def omegas(eta, **kw):
    H = 1 / eta
    oN = fe.solve(1.0, 0.3, 1.0, 0.3, H, 1.0, 0.0, **kw)
    oM = fe.solve(1.0, 0.3, 1.0, 0.3, H, 0.0, 1.0, **kw)
    return (oN['psi_deg'] + 180) % 360, (oM['psi_deg'] - gamma_deg(0.0, eta) + 90) % 360, oN['ndof']


out = []
t0 = time.time()
runs = [(eta, {'far': 0.1}) for eta in (0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)] + \
       [(0.5, {}), (0.5, {'refine': 1.5}), (0.5, {'Lfac': 14.0}), (0.5, {'fit': (2e-3, 2e-2)}), (0.5, {'fit': (1e-2, 8e-2)})]
for eta, kw in runs:
    wN, wM, nd = omegas(eta, **kw)
    out.append({'eta': eta, 'options': kw, 'omega_from_N': wN, 'omega_from_M': wM, 'fit_1992': 52.1 - 3 * eta, 'ndof': nd})
    print('eta %.2f %-28s omega N %.3f  M %.3f   1992 fit %.2f   ndof %d   %.0f s' % (eta, kw, wN, wM, 52.1 - 3 * eta, nd, time.time() - t0), flush=True)
json.dump(out, open(os.path.join(HERE, 'homogeneous.json'), 'w'), indent=1)
