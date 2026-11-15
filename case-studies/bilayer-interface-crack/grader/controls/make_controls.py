"""Two control submissions for the grader of the bilayer problem, written by the organizer, with the grades expected.

c1_own_fe          the phase angles from a finite-element solution (our first route: its code, a run script and its log
                   are in the folder); every number passes. Expected: own_computation, error none.
c2_printed_table   the closed forms right; the phase angles from the angle omega of the source paper's printed table,
                   named as the source; no computation. Expected: stated_with_source, error theory.
The numbers are read from ../../key/reference.json. usage: python3 -I make_controls.py OUT_DIR   (OUT_DIR must not exist)
"""
import json, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = os.path.join(HERE, '..', '..', 'key')
REF = json.load(open(os.path.join(KEY, 'reference.json')))
out = sys.argv[1]
if os.path.exists(out):
    sys.exit('%s exists; refusing to overwrite' % out)


def results(pa, pb, pboth):
    return {'eps': round(REF['eps'], 6),
            'caseA': {'G_E2_h_over_P2': round(REF['caseA']['G_E2_h_over_P2'], 6), 'Kabs_sqrt_h_over_P': round(REF['caseA']['Kabs_sqrt_h_over_P'], 6), 'psi_deg': round(pa, 3)},
            'caseB': {'G_E2_h3_over_M2': round(REF['caseB']['G_E2_h3_over_M2'], 6), 'Kabs_h32_over_M': round(REF['caseB']['Kabs_h32_over_M'], 6), 'psi_deg': round(pb, 3)},
            'psi_both_deg': round(pboth, 3)}


CLOSED = """## Closed forms

Units h = 1, E2 = 1. Plane strain: Ebar = E/(1 - nu^2), so Ebar1 = 81/8 and Ebar2 = 9/8; mu1 = 27/8, mu2 = 3/8,
kappa1 = kappa2 = 5/3.

- eps = (1/2 pi) ln(2/3) = -0.064532.
- Far ahead of the tip the bonded part is free of stress, so G is the energy per unit length of the two separate
  layers far behind the tip: G = P^2/(2 Ebar1 h) + 6 M^2/(Ebar1 h^3) + P^2/(2 Ebar2 H) + 6 (M + P (h + H)/2)^2/(Ebar2 H^3).
  Case A: 287/162 = 1.771605. Case B: 34/27 = 1.259259.
- Two materials: G = (c1 + c2) |K|^2/(16 cosh^2(pi eps)), c = (kappa + 1)/mu. Case A: |K| = 1.933127. Case B: 1.629801.
"""

# ---- c1: own finite-element computation, right ----
f = REF['first_route']
d = os.path.join(out, 'c1_own_fe')
os.makedirs(d)
shutil.copyfile(os.path.join(KEY, 'fe_bilayer.py'), os.path.join(d, 'fe_bilayer.py'))
open(os.path.join(d, 'run.py'), 'w').write('''"""Both load cases at three mesh densities; prints psi, J against the beam G, and |K| from the crack opening against |K| from G."""
import fe_bilayer as fe
E1, NU1, E2, NU2, H = 9.0, 1/3, 1.0, 1/3, 2.0
for refine in (1.0, 1.5, 2.0):
    a = fe.solve(E1, NU1, E2, NU2, H, -1.0, 0.0, refine=refine)      # case A: compressive force P = 1 (N = -P)
    b = fe.solve(E1, NU1, E2, NU2, H, 0.0, 1.0, refine=refine)       # case B: moment M = 1
    print('refine %.1f ndof %d  psiA %.4f  psiB %.4f  J/G %.5f %.5f  |K|cod/|K|G %.4f %.4f' % (
        refine, a['ndof'], a['psi_deg'], b['psi_deg'], a['J'] / a['G_beam'], b['J'] / b['G_beam'],
        abs(a['K']) / a['Kabs_J'], abs(b['K']) / b['Kabs_J']))
''')
log = ''.join('refine %.1f ndof %d  psiA %.4f  psiB %.4f  J/G %.5f %.5f  |K|cod/|K|G %.4f %.4f\n' % (
    r['refine'], r['ndof'], r['psiA'], r['psiB'], r['J_over_G_A'], r['J_over_G_B'], r['Kcod_over_exact_A'], r['Kcod_over_exact_B'])
    for r in f['mesh_sequence'])
open(os.path.join(d, 'run.log'), 'w').write(log)
json.dump(results(f['psiA_deg'], f['psiB_deg'], f['psi_both_deg']), open(os.path.join(d, 'results.json'), 'w'), indent=1)
open(os.path.join(d, 'method.md'), 'w').write('''# Method

''' + CLOSED + '''
## Phase angles: finite elements

No formula gives the phase angle, so I solved the plane-strain problem by finite elements (`fe_bilayer.py`, run by
`run.py`, output in `run.log`): 9-node quadrilaterals graded towards the tip, the two layers 8 (h + H) long on each side
of the tip, the stresses of the task applied as tractions at the left end. K is read from the crack opening behind the
tip, delta2 + i delta1 = (c1 + c2) K r^(i eps) sqrt(r/2 pi) / (2 (1 + 2 i eps) cosh(pi eps)), by a fit over a range of r,
and psi = arg(K h^(i eps)) with h = 1.

Finest mesh: psi_A = %.3f degrees, psi_B = %.3f degrees. Both loads, M = P h: K is linear in the loads, so I added the
two complex numbers |K_A| e^(i psi_A) + |K_B| e^(i psi_B): psi = %.3f degrees.

## Checks

1. The J integral of the finite-element solution equals the beam value of G within 2e-5 for both load cases (`run.log`,
   column J/G).
2. Three mesh densities: psi_A changes by 0.016 degree and psi_B by 0.010 degree from the coarsest to the finest
   (`run.log`).
3. The magnitude of K from the crack opening is within 0.4 %% of the magnitude from G on every mesh, and within 0.2 %% on
   the finest (`run.log`, last two columns).
''' % (f['psiA_deg'], f['psiB_deg'], f['psi_both_deg']))

# ---- c2: the printed table value, named as the source ----
t = REF['record_printed_table']
d = os.path.join(out, 'c2_printed_table')
os.makedirs(d)
json.dump(results(t['psiA_deg'], t['psiB_deg'], t['psi_both_deg']), open(os.path.join(d, 'results.json'), 'w'), indent=1)
open(os.path.join(d, 'combine.py'), 'w').write('''"""psi under both loads from the two complex stress intensity factors (per unit P, h = 1, M = P h)."""
import cmath, math
KA, KB, omega, gamma = %.6f, %.6f, %.1f, %.4f
kA = KA * cmath.exp(1j * math.radians(omega))
kB = KB * cmath.exp(1j * math.radians(omega + gamma - 90))
print(omega, omega + gamma - 90, math.degrees(cmath.phase(kA + kB)))
''' % (REF['caseA']['Kabs_sqrt_h_over_P'], REF['caseB']['Kabs_h32_over_M'], t['omega_deg'], REF['gamma_deg']))
open(os.path.join(d, 'method.md'), 'w').write('''# Method

''' + CLOSED + '''
## Phase angles

The solution of this problem is known: Z. Suo and J. W. Hutchinson, Interface crack between two elastic layers,
Int. J. Fract. 43 (1990). There K h^(i eps) = sqrt((1 - alpha)/(1 - beta^2)) (P/sqrt(2 h U) - i e^(i gamma) M/sqrt(2 h^3 V)) e^(i omega),
where gamma follows from the beam constants and omega(alpha, beta, h/H) is tabulated. Here alpha = 0.8, beta = 0.2,
h/H = 0.5, and the table gives omega = %.1f degrees. I did not compute omega myself.

- Case A: psi = omega = %.3f degrees.
- Case B: psi = omega + gamma - 90 degrees, with gamma = %.3f degrees: psi = %.3f degrees.
- Both loads, M = P h: the two complex numbers added (`combine.py`): psi = %.3f degrees.

## Checks

1. The difference psi_A - psi_B = 90 - gamma = %.3f degrees agrees with the energy: cos(psi_A - psi_B) = b/sqrt(a c)
   with a = 287/162, b = 1, c = 34/27 gives the same angle.
''' % (t['omega_deg'], t['psiA_deg'], REF['gamma_deg'], t['psiB_deg'], t['psi_both_deg'], 90 - REF['gamma_deg']))

json.dump({'c1_own_fe': {'theory_source': 'own_computation', 'error_location': 'none', 'score': '8 of 8'},
           'c2_printed_table': {'theory_source': 'stated_with_source', 'error_location': 'theory', 'score': '5 of 8'}},
          open(os.path.join(out, 'EXPECTED.json'), 'w'), indent=1)
print('controls in', out)
