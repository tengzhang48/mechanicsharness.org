"""Validation of fe_bilayer.py against results that do not come from it.

1. Homogeneous layers of equal thickness: exact K (Suo and Hutchinson 1990, Eq. 2.19), omega = arccos(sqrt(3/7)) = 49.107 deg.
2. Entries of the paper's Tables 1-4 (omega(alpha, beta, eta) in degrees, from their integral-equation solution, quoted
   to 0.1 deg): material constants are constructed for each (alpha, beta) with nu2 = 0.3 (plane strain).
3. Mesh refinement: the same case at refine = 1, 1.5, 2.
For each case, omega is obtained twice: from the axial load (psi_N + 180 deg: the source's P is minus our N) and from
the moment (psi_M - gamma + 90 deg), gamma from the beam constants (Eq. 2.8). Both use only the finite-element phase.
|K| is checked against J and against the beam energy. Writes validation.json. Run: python3 -I validate_fe.py
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fe_bilayer as fe  # noqa: E402


def materials(alpha, beta, nu2=0.3):
    k2 = 3 - 4 * nu2
    S = 2 / (beta - (1 + alpha) * (k2 - 1) / (2 * (k2 + 1)) + (1 - alpha) / 2)
    k1 = (1 - alpha) * S / 2 - 1
    Gm = (1 + alpha) * S / (2 * (k2 + 1))
    nu1 = (3 - k1) / 4
    if not (0 <= nu1 <= 0.5) or Gm <= 0:
        raise ValueError('nu1 = %.3f out of range, or a negative modulus ratio' % nu1)
    mu2 = 1 / (2 * (1 + nu2))          # E2 = 1
    mu1 = Gm * mu2
    return 2 * mu1 * (1 + nu1), nu1, 1.0, nu2


def gamma_deg(alpha, eta):
    Sig = (1 + alpha) / (1 - alpha)
    A = 1 / (1 + Sig * (4 * eta + 6 * eta ** 2 + 3 * eta ** 3))
    I = 1 / (12 * (1 + Sig * eta ** 3))
    return math.degrees(math.asin(6 * Sig * eta ** 2 * (1 + eta) * math.sqrt(A * I)))


def materials_any(alpha, beta):
    for nu2 in (0.3, 0.2, 0.4, 0.1, 0.45, 0.0, 0.49) + tuple(k / 300 for k in range(1, 150)):   # then a scan in steps of 1/300
        try:
            return materials(alpha, beta, nu2)
        except (ValueError, ZeroDivisionError):
            pass
    raise ValueError('no admissible Poisson ratios')


def omega(alpha, beta, eta, refine=1.0):
    E1, nu1, E2, nu2 = materials_any(alpha, beta)
    H = 1 / eta
    oN = fe.solve(E1, nu1, E2, nu2, H, 1.0, 0.0, refine=refine)
    oM = fe.solve(E1, nu1, E2, nu2, H, 0.0, 1.0, refine=refine)
    a, b = oN['alpha'], oN['beta']
    assert abs(a - alpha) < 1e-9 and abs(b - beta) < 1e-9
    wN = (oN['psi_deg'] + 180) % 360
    wM = (oM['psi_deg'] - gamma_deg(alpha, eta) + 90) % 360
    return {'alpha': alpha, 'beta': beta, 'eta': eta, 'refine': refine, 'E1': E1, 'nu1': nu1, 'E2': E2, 'nu2': nu2,
            'omega_from_N': wN, 'omega_from_M': wM, 'ndof': oN['ndof'],
            'J_over_Gbeam_N': oN['J'] / oN['G_beam'], 'J_over_Gbeam_M': oM['J'] / oM['G_beam'],
            'Kcod_over_KJ_N': abs(oN['K']) / oN['Kabs_J'], 'Kcod_over_KJ_M': abs(oM['K']) / oM['Kabs_J']}


CASES = [  # (alpha, beta, eta, table value or exact)
    (0.0, 0.0, 1.0, math.degrees(math.acos(math.sqrt(3 / 7)))),
    (0.0, 0.0, 0.5, 49.8), (0.0, 0.0, 0.1, 54.6),
    (0.4, 0.1, 0.5, 49.4), (0.4, 0.1, 1.0, 47.8), (-0.4, -0.1, 0.5, 51.5),
    (0.8, 0.2, 0.5, 51.4), (0.6, 0.2, 0.1, 52.4), (-0.6, -0.3, 0.1, 60.5), (0.2, -0.1, 0.5, 55.6)]

if __name__ == '__main__':
    out = []
    t0 = time.time()
    if len(sys.argv) == 4:      # one point only: python3 -I validate_fe.py ALPHA BETA ETA  (prints; validation.json is not rewritten)
        a, b, eta = (float(x) for x in sys.argv[1:])
        r = omega(a, b, eta)
        ref = [c[3] for c in CASES if (c[0], c[1], c[2]) == (a, b, eta)]
        print('alpha %+.1f beta %+.1f eta %.1f | E1/E2 %.4f nu1 %.4f nu2 %.4f | omega N %.3f  M %.3f | reference %s | J/G %.5f %.5f | Kcod/KJ %.4f %.4f | %.0f s' %
              (a, b, eta, r['E1'] / r['E2'], r['nu1'], r['nu2'], r['omega_from_N'], r['omega_from_M'], ('%.3f' % ref[0]) if ref else 'none',
               r['J_over_Gbeam_N'], r['J_over_Gbeam_M'], r['Kcod_over_KJ_N'], r['Kcod_over_KJ_M'], time.time() - t0))
        sys.exit(0)
    for (a, b, eta, ref) in CASES:
        try:
            r = omega(a, b, eta)
        except ValueError as e:
            print('skip', a, b, eta, e)
            continue
        r['reference'] = ref
        out.append(r)
        print('alpha %+.1f beta %+.1f eta %.1f | omega N %.3f  M %.3f | reference %.3f | diff %+.3f %+.3f | J/G %.5f %.5f | Kcod/KJ %.4f %.4f | %.0f s' %
              (a, b, eta, r['omega_from_N'], r['omega_from_M'], ref, r['omega_from_N'] - ref, r['omega_from_M'] - ref,
               r['J_over_Gbeam_N'], r['J_over_Gbeam_M'], r['Kcod_over_KJ_N'], r['Kcod_over_KJ_M'], time.time() - t0), flush=True)
    for refine in (1.5, 2.0):
        r = omega(0.0, 0.0, 1.0, refine=refine)
        r['reference'] = CASES[0][3]
        out.append(r)
        print('refine %.1f homogeneous eta 1: omega N %.4f M %.4f (exact %.4f), ndof %d, %.0f s' %
              (refine, r['omega_from_N'], r['omega_from_M'], CASES[0][3], r['ndof'], time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(HERE, 'validation.json'), 'w'), indent=1)
