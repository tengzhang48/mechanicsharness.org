"""Reference answers of the bilayer interface-crack benchmark, written to reference.json.

Exact parts (closed form): eps, the energy release rates (energy far behind the tip, beam theory; far ahead stress-free),
and |K| from G = (c1 + c2) |K|^2/(16 cosh^2(pi eps)).
Numerical part: the phase psi = arg(K h^(i eps)) of each load case. It was computed by two independent finite-element
routes. The key takes the second route's values (second_route/, see its REPORT.md): it is the more accurate of the two
(it meets the exact relative phase of the two loadings to 1e-5 degree and exact cases to 0.001 degree). The first
route (fe_bilayer.py, validated in validate_fe.py, check_griffith_interface.py and check_homogeneous.py) is run here at
three mesh densities as the check; its finest mesh and the differences from the key are recorded, and the script stops
if the routes differ by more than 0.1 degree. Combined quantity: psi with both loads acting, M = P h (K is linear in
the loads).
For the record (not used by the key): omega and the answers obtained from the source paper's Tables 1-4 (Suo and
Hutchinson 1990) as printed, interpolated if needed: the answer of a solver that has those tables.
Run: python3 -I reference_bilayer.py   (NumPy, SciPy; about 10 minutes)
"""
import cmath, json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fe_bilayer as fe  # noqa: E402

# Parameters of the task (changed on 13 Nov 2026; before: E1/E2 = 5, nu1 = 0.22, nu2 = 0.35, h/H = 0.3). The new set has
# Dundurs parameters alpha = 0.8, beta = 0.2 and thickness ratio 0.5: a printed entry of the source paper's Table 3
# (51.4 degrees), where our finite elements give a value about 0.9 degree higher, so that the printed entry does not
# pass the 0.5 degree check; and a bimaterial constant large enough (eps = -0.0645) for the two-material relation between
# G and |K| (2 %) and the choice of reference length (2.6 degrees) to be tested at the stated tolerances.
E1, NU1, E2, NU2 = 9.0, 1.0 / 3.0, 1.0, 1.0 / 3.0
ETA = 0.5                     # h/H
CHI_BOTH = 1.0                # the combined loading of the task: M = CHI_BOTH * P * h
SECOND_ROUTE_SET = 'S2'       # the name of these parameters in second_route/results/final.json
H_OVER_h = 1 / ETA
P0, M0 = 1.0, 1.0             # load case A: upper arm force N = -P0 (compressive); case B: moment M0 (opening)


def exact_parts():
    alpha, beta, eps, c1, c2 = fe.dundurs(E1, NU1, E2, NU2)
    GA = fe.beam_G(E1, NU1, E2, NU2, 1.0, H_OVER_h, -P0, 0.0)
    GB = fe.beam_G(E1, NU1, E2, NU2, 1.0, H_OVER_h, 0.0, M0)
    kmag = lambda G: math.sqrt(16 * math.cosh(math.pi * eps) ** 2 * G / (c1 + c2))
    return alpha, beta, eps, c1, c2, GA, GB, kmag(GA), kmag(GB)


def table_omega(alpha, beta, eta):
    """Trilinear interpolation of Suo and Hutchinson (1990), Tables 1-4 (degrees), as printed."""
    A = [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]
    B = [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
    # only the cells around (alpha, beta) = (0.6..0.8, 0.1..0.2) are needed here: rows beta = 0.1, 0.2; columns 0.6, 0.8
    T = {0.0: {0.1: {0.6: 56.1, 0.8: 62.0}, 0.2: {0.6: 52.1, 0.8: 57.0}},
         0.1: {0.1: {0.6: 56.5, 0.8: 60.9}, 0.2: {0.6: 52.4, 0.8: 57.7}},
         0.5: {0.1: {0.6: 51.8, 0.8: 55.2}, 0.2: {0.6: 47.8, 0.8: 51.4}},
         1.0: {0.1: {0.6: 49.8, 0.8: 52.4}, 0.2: {0.6: 45.1, 0.8: 47.6}}}
    def bil(t):
        ta = (alpha - 0.6) / 0.2; tb = (beta - 0.1) / 0.1
        return ((1 - ta) * (1 - tb) * t[0.1][0.6] + ta * (1 - tb) * t[0.1][0.8] + (1 - ta) * tb * t[0.2][0.6] + ta * tb * t[0.2][0.8])
    etas = sorted(T)
    for e0, e1 in zip(etas[:-1], etas[1:]):
        if e0 <= eta <= e1:
            w = (eta - e0) / (e1 - e0)
            return (1 - w) * bil(T[e0]) + w * bil(T[e1])


def gamma_rad(alpha, eta):
    Sig = (1 + alpha) / (1 - alpha)
    A = 1 / (1 + Sig * (4 * eta + 6 * eta ** 2 + 3 * eta ** 3))
    I = 1 / (12 * (1 + Sig * eta ** 3))
    return math.asin(6 * Sig * eta ** 2 * (1 + eta) * math.sqrt(A * I)), A, I


if __name__ == '__main__':
    t0 = time.time()
    alpha, beta, eps, c1, c2, GA, GB, KA_abs, KB_abs = exact_parts()
    runs = []
    for refine in (1.0, 1.5, 2.0):
        oA = fe.solve(E1, NU1, E2, NU2, H_OVER_h, -P0, 0.0, refine=refine)
        oB = fe.solve(E1, NU1, E2, NU2, H_OVER_h, 0.0, M0, refine=refine)
        runs.append({'refine': refine, 'ndof': oA['ndof'], 'psiA': oA['psi_deg'], 'psiB': oB['psi_deg'],
                     'J_over_G_A': oA['J'] / GA, 'J_over_G_B': oB['J'] / GB,
                     'Kcod_over_exact_A': abs(oA['K']) / KA_abs, 'Kcod_over_exact_B': abs(oB['K']) / KB_abs})
        print('refine %.1f ndof %d  psiA %.4f  psiB %.4f  J/G %.5f %.5f  |K|cod/exact %.4f %.4f  %.0f s' %
              (refine, oA['ndof'], oA['psi_deg'], oB['psi_deg'], runs[-1]['J_over_G_A'], runs[-1]['J_over_G_B'],
               runs[-1]['Kcod_over_exact_A'], runs[-1]['Kcod_over_exact_B'], time.time() - t0), flush=True)
    both = lambda a, b: math.degrees(cmath.phase(KA_abs * cmath.exp(1j * math.radians(a)) + CHI_BOTH * KB_abs * cmath.exp(1j * math.radians(b))))
    psiA1, psiB1 = runs[-1]['psiA'], runs[-1]['psiB']       # first route, finest mesh
    # the key: the second route (see the docstring)
    r2 = json.load(open(os.path.join(HERE, 'second_route', 'results', 'final.json')))[SECOND_ROUTE_SET]
    assert abs(r2['A']['eps'] - eps) < 1e-9 and abs(r2['A']['absK_G'] - KA_abs) < 1e-6 and abs(r2['B']['absK_G'] - KB_abs) < 1e-6, \
        'the second route was run for other parameters'
    psiA, psiB = r2['A']['psi'], r2['B']['psi']
    psi_both = both(psiA, psiB)     # both loads with M = CHI_BOTH P h: K h^(i eps) = K_A + CHI_BOTH K_B (K is linear in the loads)
    gam, Ab, Ib = gamma_rad(alpha, ETA)
    omega_fe = (psiA - 0.0) % 360                          # case A has P = +P0 in the source's convention: psi_A = omega
    omega_fe_B = (psiB - math.degrees(gam) + 90) % 360
    om_tab = table_omega(alpha, beta, ETA)
    kAt = KA_abs * cmath.exp(1j * math.radians(om_tab)); kBt = KB_abs * cmath.exp(1j * math.radians(om_tab + math.degrees(gam) - 90))
    ref = {'materials': {'E1_over_E2': E1 / E2, 'nu1': NU1, 'nu2': NU2, 'h_over_H': ETA},
           'eps': eps, 'alpha': alpha, 'beta': beta,
           'caseA': {'G_E2_h_over_P2': GA, 'Kabs_sqrt_h_over_P': KA_abs, 'psi_deg': psiA},
           'caseB': {'G_E2_h3_over_M2': GB, 'Kabs_h32_over_M': KB_abs, 'psi_deg': psiB},
           'phase_angles_from': 'second_route, set %s' % SECOND_ROUTE_SET,
           'relative_phase_deg': {'exact_from_energy': 90 - math.degrees(gam), 'key': psiA - psiB, 'first_route': psiA1 - psiB1},
           'M_over_Ph_both': CHI_BOTH, 'psi_both_deg': psi_both,
           'omega_from_A_deg': omega_fe, 'omega_from_B_deg': omega_fe_B, 'gamma_deg': math.degrees(gam),
           'record_printed_table': {'omega_deg': om_tab, 'psiA_deg': om_tab, 'psiB_deg': om_tab + math.degrees(gam) - 90,
                                    'psi_both_deg': math.degrees(cmath.phase(kAt + CHI_BOTH * kBt))}}
    # the first route as the check: its finest mesh, the differences from the key (first route minus key), its mesh sequence
    ref['first_route'] = {'psiA_deg': psiA1, 'psiB_deg': psiB1, 'psi_both_deg': both(psiA1, psiB1),
                          'diff_psiA_deg': psiA1 - psiA, 'diff_psiB_deg': psiB1 - psiB, 'diff_psi_both_deg': both(psiA1, psiB1) - psi_both,
                          'mesh_sequence': runs}
    assert max(abs(ref['first_route'][k]) for k in ('diff_psiA_deg', 'diff_psiB_deg', 'diff_psi_both_deg')) < 0.1, \
        'the two routes differ by more than 0.1 degree: %r' % ref['first_route']
    json.dump(ref, open(os.path.join(HERE, 'reference.json'), 'w'), indent=1)
    print(json.dumps(ref, indent=1))
