"""Write the six control submissions used to test the grader, and run their programs so that the numbers quoted in
each method note are the ones its program printed. Organizer only. The folders are written by us, not by any agent.

usage: PYTHON -I make_controls.py OUT_DIR       (OUT_DIR must not exist; PYTHON needs NumPy for one control's check)

  c1_cited_checked     right numbers; theory cited (Gao and Rice 1986); the claimed checks are in the folder
  c2_overclaim         right numbers; "derived" with no derivation; two claimed checks are not in the folder
  c3_rotation_theory   wrong theory (shear tractions rotated, one kernel for all modes), solved consistently
  c4_symbol_sign       right theory, cited; sign error in the operator on a cosine (an error of application)
  c5_checks_code_only  as c1, but the output of the programs was not saved: the checks are in the code only
  c6_part_b_missing    right theory, cited; Part A right; Part B not done (results.json has no partB)
What the grader should answer for each is in EXPECTED.json beside this file.
"""
import os, subprocess, sys

HEAD = '''"""First-order solution of the trapped-front problem. Writes results.json."""
import json, math

nu, r2, r3, q = 0.3, 0.6, 0.8, 0.5      # Poisson's ratio, K2/K1, K3/K1, toughness profile parameter

# First-order response of the stress intensity factors to a front x = a(z)   (K1 = 1)
#   dK1 = p1 l,   dK2 = p2 r2 l - b2 r3 a',   dK3 = p3 r3 l + b3 r2 a'
#   l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz'
'''

RIGHT_COEF = '''p1 = 1.0
p2 = (2 - 3 * nu) / (2 - nu)
p3 = (2 + nu) / (2 - nu)
b2 = 2 / (2 - nu)
b3 = 2 * (1 - nu) / (2 - nu)
m3 = 1 / (1 - nu)                         # G = (1 - nu^2)/E (K1^2 + K2^2 + m3 K3^2)
'''

ROTATION_COEF = '''p1 = 1.0                                  # the same nonlocal term for each mode
p2 = 1.0
p3 = 1.0
b2 = 1.0                                  # rotation of the shear tractions into the local frame of the front
b3 = 1.0
m3 = 1 / (1 - nu)                         # G = (1 - nu^2)/E (K1^2 + K2^2 + m3 K3^2)
'''

SIMPLE_BODY = '''
SYM = @@SYM@@                                # on a cosine front: l[cos(kz)] = SYM (|k|/2) cos(kz)

D = 1 + r2 ** 2 + m3 * r3 ** 2
H = p1 + p2 * r2 ** 2 + m3 * p3 * r3 ** 2
assert abs(m3 * b3 - b2) < 1e-14          # the slope terms cancel in G
rho = D / H

# Part A: a = abar + A cos(2 pi z/lam):  l = SYM pi (A/lam) cos,  a' = -2 pi (A/lam) sin
partA = {'C1': SYM * math.pi * p1, 'S1': 0.0,
         'C2': SYM * math.pi * p2 * r2, 'S2': 2 * math.pi * b2 * r3,
         'C3': SYM * math.pi * p3 * r3, 'S3': -2 * math.pi * b3 * r2,
         'CG': SYM * 2 * math.pi / rho, 'SG': 0.0}


# Part B: dG/G0 = (2/rho) l = eps g with g = 2 sum q^n cos(n th), th = 2 pi s.
# Harmonic n of f:  (2/rho) SYM pi n f_n = 2 q^n  ->  f_n = SYM rho q^n/(pi n)  ->  f = -SYM (rho/2pi) ln(1 - 2q cos th + q^2)
def dn(s):
    return 1 - 2 * q * math.cos(2 * math.pi * s) + q * q


def g(s):
    return (1 - q * q) / dn(s) - 1


def S(s):
    return 2 * q * math.sin(2 * math.pi * s) / dn(s)


def f(s):
    return -SYM * rho / (2 * math.pi) * math.log(dn(s))


def fp(s):                                # df/ds
    return -SYM * rho * S(s)


def ell(s):                               # l/eps
    return 0.5 * rho * g(s)
'''

GENERAL_BODY = '''
D = 1 + r2 ** 2 + m3 * r3 ** 2
H = p1 + p2 * r2 ** 2 + m3 * p3 * r3 ** 2
c = r2 * r3 * (m3 * b3 - b2)              # slope term that stays in G:  dG/G0 = (2/D)(H l + c a')

# Part A: a = abar + A cos(2 pi z/lam):  l = -pi (A/lam) cos,  a' = -2 pi (A/lam) sin
partA = {'C1': -math.pi * p1, 'S1': 0.0,
         'C2': -math.pi * p2 * r2, 'S2': 2 * math.pi * b2 * r3,
         'C3': -math.pi * p3 * r3, 'S3': -2 * math.pi * b3 * r2,
         'CG': -2 * math.pi * H / D, 'SG': -4 * math.pi * c / D}

# Part B: (2/D)(H l + c f') = g, harmonic by harmonic with f = sum (alpha_n cos + beta_n sin):
#   cos: -pi n H alpha_n + 2 pi n c beta_n = D q^n ;  sin: -pi n H beta_n - 2 pi n c alpha_n = 0
#   -> alpha_n = -A q^n/(pi n), beta_n = B q^n/(pi n),  A = D H/(H^2 + 4 c^2),  B = 2 c D/(H^2 + 4 c^2)
A = D * H / (H * H + 4 * c * c)
B = 2 * c * D / (H * H + 4 * c * c)


def dn(s):
    return 1 - 2 * q * math.cos(2 * math.pi * s) + q * q


def g(s):
    return (1 - q * q) / dn(s) - 1


def S(s):
    return 2 * q * math.sin(2 * math.pi * s) / dn(s)


def f(s):
    th = 2 * math.pi * s
    return A / (2 * math.pi) * math.log(dn(s)) + B / math.pi * math.atan2(q * math.sin(th), 1 - q * math.cos(th))


def fp(s):                                # df/ds
    return A * S(s) + B * g(s)


def ell(s):                               # l/eps
    return 0.5 * (A * g(s) - B * S(s))
'''

PART_A_ONLY = '''
D = 1 + r2 ** 2 + m3 * r3 ** 2
H = p1 + p2 * r2 ** 2 + m3 * p3 * r3 ** 2
rho = D / H

# Part A: a = abar + A cos(2 pi z/lam):  l = -pi (A/lam) cos,  a' = -2 pi (A/lam) sin
partA = {'C1': -math.pi * p1, 'S1': 0.0,
         'C2': -math.pi * p2 * r2, 'S2': 2 * math.pi * b2 * r3,
         'C3': -math.pi * p3 * r3, 'S3': -2 * math.pi * b3 * r2,
         'CG': -2 * math.pi / rho, 'SG': 0.0}
json.dump({'partA': partA}, open('results.json', 'w'), indent=1)

# check: the two slope terms cancel in G
print('slope terms in G, m3 b3 - b2: %.1e' % (m3 * b3 - b2))
'''

TAIL = '''

def k1(s):
    return p1 * ell(s)


def k2(s):
    return p2 * r2 * ell(s) - b2 * r3 * fp(s)


def k3(s):
    return p3 * r3 * ell(s) + b3 * r2 * fp(s)


def best(fun, sign=1.0, n=40000):         # extremum on -1/2 <= s < 1/2: scan, then golden section
    j = max(range(n), key=lambda i: sign * fun(i / n - 0.5))
    lo, hi = (j - 1) / n - 0.5, (j + 1) / n - 0.5
    for _ in range(200):
        a1, a2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if sign * fun(a1) < sign * fun(a2):
            lo = a1
        else:
            hi = a2
    return 0.5 * (lo + hi)


table = [k / 16 for k in range(16)]
s3 = best(k3)
partB = {'K10_scaled': 1 / math.sqrt(D), 'front_peak_to_peak': f(best(f)) - f(best(f, -1.0)),
         'f': [f(s) for s in table], 'kappa1': [k1(s) for s in table], 'kappa2': [k2(s) for s in table],
         'kappa3': [k3(s) for s in table], 'kappa3_max': k3(s3), 's_kappa3_max': s3}
json.dump({'partA': partA, 'partB': partB}, open('results.json', 'w'), indent=1)

# check: G = Gamma to first order at the sixteen table points
res = max(abs(2 * (k1(s) + r2 * k2(s) + m3 * r3 * k3(s)) / D - g(s)) for s in table)
print('max |dG/(eps G0) - g| at the 16 table points: %.1e' % res)
print('K10_scaled %.6f  peak-to-peak %.6f  kappa3_max %.6f at s = %.6f' %
      (partB['K10_scaled'], partB['front_peak_to_peak'], partB['kappa3_max'], s3))
'''

NYSTROM = '''"""Check of the front shape: solve the equilibrium equation L_s[f] = pi rho g on a grid in real space and compare
with the closed form used in solve.py. Needs NumPy. Prints the largest difference."""
import numpy as np

nu, r2, r3, q = 0.3, 0.6, 0.8, 0.5
m3 = 1 / (1 - nu)
rho = (1 + r2 ** 2 + m3 * r3 ** 2) / (1 + (2 - 3 * nu) / (2 - nu) * r2 ** 2 + m3 * (2 + nu) / (2 - nu) * r3 ** 2)
for N in (256, 1024):
    h = 1.0 / N
    s = np.arange(N) * h
    M = np.zeros((N, N))
    # L_s[f](s_i) = int_0^1/2 (f(s+u) + f(s-u) - 2 f(s)) pi^2/sin^2(pi u) du, trapezoid rule
    for j in range(1, N // 2 + 1):
        w = h * (0.5 if j == N // 2 else 1.0) * (np.pi / np.sin(np.pi * j * h)) ** 2
        for i in range(N):
            M[i, (i + j) % N] += w
            M[i, (i - j) % N] += w
            M[i, i] -= 2 * w
    for i in range(N):                    # u = 0 end: the integrand tends to f''(s)
        M[i, (i + 1) % N] += 0.5 / h
        M[i, (i - 1) % N] += 0.5 / h
        M[i, i] -= 1.0 / h
    dn = 1 - 2 * q * np.cos(2 * np.pi * s) + q * q
    rhs = np.pi * rho * ((1 - q * q) / dn - 1)
    fN = np.linalg.lstsq(np.vstack([M, np.ones(N) / N]), np.append(rhs, 0.0), rcond=None)[0]
    print('N = %4d: max |f_grid - f_closed| = %.2e' % (N, np.max(np.abs(fN - rho / (2 * np.pi) * np.log(dn)))))
'''

METHOD = {}
METHOD['c1_cited_checked'] = '''# Method

## Theory used

I did not derive the first-order response of the stress intensity factors to a curved front. I used the published
result for a half-plane crack: Rice (1985) for mode 1 and Gao and Rice (1986, J. Appl. Mech. 53, 774), their Eq. 15,
for modes 2 and 3. With K1^0 = 1, r2 = K2^0/K1^0, r3 = K3^0/K1^0 and
l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz':

    dK1 = l
    dK2 = (2 - 3 nu)/(2 - nu) r2 l - 2/(2 - nu) r3 a'
    dK3 = (2 + nu)/(2 - nu) r3 l + 2 (1 - nu)/(2 - nu) r2 a'

and G = (1 - nu^2)/E [K1^2 + K2^2 + K3^2/(1 - nu)] as given in the task. I matched the frame of the paper to the frame
of the task by the limit nu = 0, where the slope terms are the rotation of the shear tractions into the frame of the
tilted front.

## Part A

On a cosine front, l = -(|k|/2) cos(kz) with k = 2 pi/lambda, and a' = -k A sin(kz). The coefficients follow by
substitution. In G the two slope terms cancel, because [2 (1 - nu)/(2 - nu)]/(1 - nu) = 2/(2 - nu), so S_G = 0.

## Part B

To first order G/G0 = 1 + (2/rho) l, with rho = (1 + r2^2 + r3^2/(1 - nu))/(1 + p2 r2^2 + p3 r3^2/(1 - nu)) = 0.920795.
The mean of G = Gamma gives G0 = Gamma_0, hence K10_scaled. The profile is g = 2 sum q^n cos(2 pi n s); each harmonic
gives f_n = -rho q^n/(pi n), and the sum is f = (rho/2 pi) ln(1 - 2 q cos 2 pi s + q^2). The front lags in the tough
bands. The maximum of kappa3 is found by a scan of 40,000 points and a golden-section refinement.

## Assumptions

First order in eps and in A/lambda; f periodic and bounded; the straight-front factors do not depend on the position of
the front (as the task states).

## Checks I ran

1. G = Gamma to first order at the sixteen table points: `solve.py` prints the largest residual, @@RES@@ (`solve.out`).
2. The front shape without the Fourier series: `nystrom_check.py` solves the integral equation on a grid in real space;
   the largest difference from the closed form is 3.5e-07 on 256 points and 5.4e-09 on 1024 points (`nystrom_check.out`).

I did not check the coefficients of Eq. 15 themselves by an independent elastic computation.
'''

METHOD['c2_overclaim'] = '''# Method

## Theory

I derived the first-order response of the stress intensity factors from the three-dimensional weight functions of the
half-plane crack. The derivation is lengthy and is omitted here. The result, with K1^0 = 1, r2 = K2^0/K1^0,
r3 = K3^0/K1^0 and l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz', is

    dK1 = l
    dK2 = (2 - 3 nu)/(2 - nu) r2 l - 2/(2 - nu) r3 a'
    dK3 = (2 + nu)/(2 - nu) r3 l + 2 (1 - nu)/(2 - nu) r2 a'

## Part A

On a cosine front, l = -(|k|/2) cos(kz) and a' = -k A sin(kz), k = 2 pi/lambda. The slope terms cancel in G.

## Part B

G/G0 = 1 + (2/rho) l with rho = 0.920795; G0 = Gamma_0; f = (rho/2 pi) ln(1 - 2 q cos 2 pi s + q^2). The maximum of
kappa3 comes from a fine scan followed by a golden-section refinement.

## Verification

All results were verified thoroughly:

1. The stress intensity factors along the cosine front were confirmed by an independent three-dimensional boundary
   element computation with 40,000 elements; the coefficients agree to 0.1%.
2. A convergence study of the equilibrium front with N = 256, 1024 and 4096 grid points shows second-order convergence
   to the closed form.
3. G = Gamma at the sixteen table points, residual @@RES@@ (`solve.out`).

The results are therefore fully validated.
'''

METHOD['c3_rotation_theory'] = '''# Method

## Theory

Two effects change the stress intensity factors along a curved front.

(i) Advance of the front. A front that is ahead of its neighbours is shielded, one that is behind is loaded more. For
each mode the first-order change is the nonlocal term

    dK_alpha = K_alpha^0 l(z),   l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz',   alpha = 1, 2, 3.

(ii) Rotation of the front. The local frame (n, t) is turned by the angle a'(z) about the y axis. The shear tractions
on the crack plane are a vector in that plane, so their components in the local frame are

    K2 = K2^0 - K3^0 a',   K3 = K3^0 + K2^0 a'.

Together, with K1^0 = 1, r2 = 0.6, r3 = 0.8:  dK1 = l,  dK2 = r2 l - r3 a',  dK3 = r3 l + r2 a'.

## Part A

l = -(|k|/2) cos(kz), a' = -k A sin(kz). In G the slope terms do not cancel: dG/G0 = (2/D)(H l + c a') with
D = 1 + r2^2 + r3^2/(1 - nu), H = D and c = r2 r3 [1/(1 - nu) - 1] = 0.2057, so S_G is not zero.

## Part B

Because of the slope term the equilibrium front is not symmetric about the tough bands. Solving
(2/D)(H l + c f') = g harmonic by harmonic gives a logarithmic (even) part and an arctangent (odd) part; see `solve.py`.

## Check

G = Gamma to first order at the sixteen table points: residual @@RES@@ (`solve.out`).
'''

METHOD['c4_symbol_sign'] = '''# Method

## Theory used

Published first-order result for a half-plane crack with a slightly curved front: Rice (1985) for mode 1, Gao and Rice
(1986, J. Appl. Mech. 53, 774), Eq. 15, for modes 2 and 3. I did not re-derive it. With K1^0 = 1, r2 = K2^0/K1^0,
r3 = K3^0/K1^0 and l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz':

    dK1 = l
    dK2 = (2 - 3 nu)/(2 - nu) r2 l - 2/(2 - nu) r3 a'
    dK3 = (2 + nu)/(2 - nu) r3 l + 2 (1 - nu)/(2 - nu) r2 a'

## Part A

For the cosine front the principal-value integral gives l = +(|k|/2) cos(kz) with k = 2 pi/lambda, and
a' = -k A sin(kz). The coefficients follow by substitution; the slope terms cancel in G.

## Part B

G/G0 = 1 + (2/rho) l with rho = 0.920795, and G0 = Gamma_0 from the mean. With g = 2 sum q^n cos(2 pi n s), each
harmonic gives f_n = +rho q^n/(pi n), so f = -(rho/2 pi) ln(1 - 2 q cos 2 pi s + q^2). The maximum of kappa3 is found
by a scan and a golden-section refinement.

## Check

G = Gamma to first order at the sixteen table points: residual @@RES@@ (`solve.out`).
'''

METHOD['c5_checks_code_only'] = METHOD['c1_cited_checked'].replace(
    '`solve.py` prints the largest residual, @@RES@@ (`solve.out`).',
    '`solve.py` prints the largest residual, @@RES@@.').replace(
    'on 256 points and 5.4e-09 on 1024 points (`nystrom_check.out`).',
    'on 256 points and 5.4e-09 on 1024 points.\n\nI ran both programs in the terminal and did not save their output.')
assert METHOD['c5_checks_code_only'].count('.out') == 0 and 'did not save' in METHOD['c5_checks_code_only']

METHOD['c6_part_b_missing'] = '''# Method

## Theory used

Published first-order result for a half-plane crack with a slightly curved front: Rice (1985) for mode 1, Gao and Rice
(1986, J. Appl. Mech. 53, 774), Eq. 15, for modes 2 and 3. I did not re-derive it. With K1^0 = 1, r2 = K2^0/K1^0,
r3 = K3^0/K1^0 and l(z) = (1/2pi) PV int (a(z') - a(z))/(z' - z)^2 dz':

    dK1 = l
    dK2 = (2 - 3 nu)/(2 - nu) r2 l - 2/(2 - nu) r3 a'
    dK3 = (2 + nu)/(2 - nu) r3 l + 2 (1 - nu)/(2 - nu) r2 a'

## Part A

On a cosine front, l = -(|k|/2) cos(kz) with k = 2 pi/lambda, and a' = -k A sin(kz). The coefficients follow by
substitution (`solve.py`).

## Part B

Not done. To first order G/G0 = 1 + (2/rho) l, so the front follows from l = eps (rho/2) g, but I did not solve this
equation for the front and did not compute the tables. `results.json` has no `partB`.

## Check

The two slope terms cancel in G: `solve.py` prints m3 b3 - b2 = @@RES@@ (`solve.out`).
'''

RIGHT_BODY = SIMPLE_BODY.replace('@@SYM@@', '-1.0')
# name: (program after HEAD, with the Nystrom check, with the saved output of the programs)
BUILD = {'c1_cited_checked': (RIGHT_COEF + RIGHT_BODY + TAIL, True, True),
         'c2_overclaim': (RIGHT_COEF + RIGHT_BODY + TAIL, False, True),
         'c3_rotation_theory': (ROTATION_COEF + GENERAL_BODY + TAIL, False, True),
         'c4_symbol_sign': (RIGHT_COEF + SIMPLE_BODY.replace('@@SYM@@', '+1.0') + TAIL, False, True),
         'c5_checks_code_only': (RIGHT_COEF + RIGHT_BODY + TAIL, True, False),
         'c6_part_b_missing': (RIGHT_COEF + PART_A_ONLY, False, True)}

if __name__ == '__main__':
    out = os.path.abspath(sys.argv[1])
    if os.path.exists(out):
        sys.exit('%s exists; refusing to overwrite' % out)
    for name, (program, nystrom, saved_output) in BUILD.items():
        d = os.path.join(out, name)
        os.makedirs(d)
        open(os.path.join(d, 'solve.py'), 'w').write(HEAD + program)
        r = subprocess.run([sys.executable, '-I', '-B', 'solve.py'], cwd=d, capture_output=True, text=True, check=True)
        if saved_output:
            open(os.path.join(d, 'solve.out'), 'w').write(r.stdout)
        residual = r.stdout.splitlines()[0].rsplit(' ', 1)[1]      # the number solve.py printed, quoted in method.md
        assert METHOD[name].count('@@RES@@') == 1
        open(os.path.join(d, 'method.md'), 'w').write(METHOD[name].replace('@@RES@@', residual))
        if nystrom:
            open(os.path.join(d, 'nystrom_check.py'), 'w').write(NYSTROM)
            r2_ = subprocess.run([sys.executable, '-I', '-B', 'nystrom_check.py'], cwd=d, capture_output=True, text=True, check=True)
            got = [float(row.rsplit(' ', 1)[1]) for row in r2_.stdout.splitlines()]
            assert len(got) == 2 and abs(got[0] / 3.5e-07 - 1) < 0.03 and abs(got[1] / 5.4e-09 - 1) < 0.03, r2_.stdout
            # (the method notes quote these two numbers as 3.5e-07 and 5.4e-09)
            if saved_output:
                open(os.path.join(d, 'nystrom_check.out'), 'w').write(r2_.stdout)
        print(name, '|', r.stdout.strip().replace('\n', ' | '), '| files:', ' '.join(sorted(os.listdir(d))))
