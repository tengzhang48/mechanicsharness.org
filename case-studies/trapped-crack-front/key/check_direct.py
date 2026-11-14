"""Second route to the reference of the trapped-front benchmark. It shares no code with reference.py and solves the
equilibrium equation in real space instead of summing a Fourier series.

What it does not test: the elastic theory. Both routes take the response coefficients from Gao and Rice's Eq. 15
(checked against the printed page and against Lazarus 2011, Eqs. 15-17; see the README). A misread coefficient would
pass here.

1. The operator L[a](z) = PV int (a(z') - a(z))/(z' - z)^2 dz' on cos(kz): adaptive quadrature gives -pi |k| cos(kz)
   (the symbol behind Part A).
2. The equilibrium equation of Part B in units lambda = eps = 1, L_s[f] = pi rho g(s), is solved by a Nystrom method in
   real space: the periodic kernel pi^2/sin^2(pi u) (the sum of 1/(u + m)^2 over the images), paired symmetrically,
   trapezoid rule on N points, the u = 0 value from a finite-difference second derivative, zero mean imposed. The
   solution is compared with the closed-form f of reference.py.
3. The closed-form f is put into the operator by adaptive quadrature at points off the grid (residual).
4. The Nystrom front is put into the operator again by a different rule, the symbol of step 1 applied by FFT. The
   residual L_s[f_N] - pi rho g measures the real-space solution with an operator it was not solved with.
5. K1, K2, K3 from Eq. 15 with the Nystrom front (L_s[f_N] from step 4, slope by FFT differentiation) against the
   reference tables; the peak-to-peak of the Nystrom front and the maximum of its kappa3 (parabola through the three
   grid values at the largest) against the reference.
6. Irwin's full G from K0 + eps K's of step 5 at eps = 0.01 and 0.001: G/G0 - (1 + eps g) is of second order in eps
   (the first-order terms, slope terms included, are right) and falls a hundredfold.
7. The maximum of kappa3 from its stationarity condition in closed form against the reference's numerical search.
8. K10_scaled written out again.
Run: python3 -I check_direct.py   (NumPy, SciPy). Exit 1 if a check fails.
"""
import json, math, os, sys
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(HERE, 'reference.json')))
nu, r2, r3, q = ref['nu'], ref['K20_over_K10'], ref['K30_over_K10'], ref['q']
fails = 0


def report(name, ok, detail):
    global fails
    fails += 0 if ok else 1
    print(('ok   ' if ok else 'FAIL ') + name + ': ' + detail)


# rho from the definitions, written out again here
G0n = 1 + r2 ** 2 + r3 ** 2 / (1 - nu)
Fn = 1 + (2 - 3 * nu) / (2 - nu) * r2 ** 2 + (2 + nu) / ((2 - nu) * (1 - nu)) * r3 ** 2
rho = G0n / Fn
report('rho', abs(rho - ref['derived']['rho_Gamma0_over_F']) < 1e-12, '%.12f' % rho)


def g(s):
    return (1 - q * q) / (1 - 2 * q * np.cos(2 * np.pi * s) + q * q) - 1


# 1. symbol of L on cos(kz), infinite line, symmetric pairing
for k in (1.0, 2.0, 5.0):
    z = 0.37
    def integrand(u):
        return (math.cos(k * (z + u)) + math.cos(k * (z - u)) - 2 * math.cos(k * z)) / (u * u)
    val = sum(quad(integrand, a, b, limit=20000)[0] for a, b in ((0, 1), (1, 50), (50, 2000))) + \
        (-2 * math.cos(k * z)) * (1 / 2000.0)   # tail: the cosines average out, -2cos(kz)/u^2 integrates to -2cos(kz)/2000
    report('L[cos kz] = -pi k cos kz, k = %g' % k, abs(val + math.pi * k * math.cos(k * z)) < 2e-3 * math.pi * k,
           'quadrature %.6f, symbol %.6f' % (val, -math.pi * k * math.cos(k * z)))

# 2. Nystrom solution of L_s[f] = pi rho g on one period
for N in (256, 1024):
    h = 1.0 / N
    s = np.arange(N) * h
    M = np.zeros((N, N))
    # L_s[f](s_i) = int_0^1/2 (f(s+u) + f(s-u) - 2 f(s)) pi^2/sin^2(pi u) du, trapezoid on u = j h, j = 0..N/2
    for j in range(1, N // 2 + 1):
        wgt = h * (0.5 if j == N // 2 else 1.0) * (np.pi / np.sin(np.pi * j * h)) ** 2
        for i in range(N):
            M[i, (i + j) % N] += wgt
            M[i, (i - j) % N] += wgt
            M[i, i] -= 2 * wgt
    # u = 0 end of the trapezoid: integrand -> f''(s), weight h/2; f'' by central differences
    for i in range(N):
        c = 0.5 * h / (h * h)
        M[i, (i + 1) % N] += c
        M[i, (i - 1) % N] += c
        M[i, i] -= 2 * c
    rhs = np.pi * rho * g(s)
    A = np.vstack([M, np.ones(N) / N])
    f = np.linalg.lstsq(A, np.append(rhs, 0.0), rcond=None)[0]
    f_closed = (rho / (2 * np.pi)) * np.log(1 - 2 * q * np.cos(2 * np.pi * s) + q * q)
    err = np.max(np.abs(f - f_closed))
    report('Nystrom front, N = %d, against the closed form' % N, err < (2e-4 if N == 256 else 2e-5),
           'max |f - f_closed| = %.2e (front amplitude %.3f)' % (err, f_closed.max() - f_closed.min()))
    if N == 1024:
        fN, sN = f, s

# 3. residual of the closed form by adaptive quadrature, at points off the grid
def f_cl(x):
    return (rho / (2 * math.pi)) * math.log(1 - 2 * q * math.cos(2 * math.pi * x) + q * q)
worst = 0.0
for x in (0.013, 0.071, 0.2, 0.333, 0.49):
    def integ(u):
        return (f_cl(x + u) + f_cl(x - u) - 2 * f_cl(x)) * (math.pi / math.sin(math.pi * u)) ** 2
    val = quad(integ, 0, 0.5, limit=400, epsabs=1e-12)[0]
    worst = max(worst, abs(val - math.pi * rho * float(g(x))))
report('closed-form front satisfies the equilibrium equation (quadrature)', worst < 1e-7, 'worst residual %.2e' % worst)

# 4. the Nystrom front in the operator by another rule: the symbol -pi |k| applied by FFT (k = 2 pi n per unit s)
h = 1.0 / len(fN)
kf = np.fft.rfftfreq(len(fN), d=h) * 2 * np.pi
fhat = np.fft.rfft(fN)
Lf = np.fft.irfft(-np.pi * kf * fhat, n=len(fN))
res = np.max(np.abs(Lf - np.pi * rho * g(sN)))
report('Nystrom front in the operator by its Fourier symbol', res < 1e-3 * np.pi * rho,
       'max |L_s[f_N] - pi rho g| = %.2e (right-hand side up to %.2f)' % (res, np.pi * rho * g(sN).max()))

# 5. K's from Eq. 15 with the Nystrom front
slope = np.fft.irfft(1j * kf * fhat, n=len(fN))                   # df/ds = a'(z)/eps
k1 = (1 / (2 * np.pi)) * Lf
k2 = -2 / (2 - nu) * r3 * slope + (2 - 3 * nu) / (2 - nu) * r2 * (1 / (2 * np.pi)) * Lf
k3 = 2 * (1 - nu) / (2 - nu) * r2 * slope + (2 + nu) / (2 - nu) * r3 * (1 / (2 * np.pi)) * Lf
idx = [int(round(t * len(fN))) % len(fN) for t in ref['partB']['s_table']]
for name, arr in (('kappa1', k1), ('kappa2', k2), ('kappa3', k3), ('f', fN)):
    e = max(abs(arr[i] - v) for i, v in zip(idx, ref['partB'][name]))
    report('%s table: direct route against reference' % name, e < 1e-4, 'max difference %.2e' % e)
p2p = fN.max() - fN.min()
report('front peak-to-peak: direct route against reference', abs(p2p - ref['partB']['front_peak_to_peak']) < 1e-6,
       '%.9f against %.9f' % (p2p, ref['partB']['front_peak_to_peak']))
j = int(np.argmax(k3))
ym, y0, yp = k3[j - 1], k3[j], k3[(j + 1) % len(k3)]
dj = 0.5 * (ym - yp) / (ym - 2 * y0 + yp)                         # vertex of the parabola through three grid values
s_dir, k_dir = ((sN[j] + dj * h + 0.5) % 1.0) - 0.5, y0 - 0.25 * (ym - yp) * dj
report('maximum of kappa3: direct route against reference',
       abs(k_dir - ref['partB']['kappa3_max']) < 1e-4 and abs(s_dir - ref['partB']['s_kappa3_max']) < 1e-4,
       'kappa3_max %.7f at s = %.7f against %.7f at %.7f' % (k_dir, s_dir, ref['partB']['kappa3_max'], ref['partB']['s_kappa3_max']))

# 6. Irwin's full G at finite eps: G/G0 - (1 + eps g) must be of second order in eps
errs = []
for eps in (1e-2, 1e-3):
    Gfull = ((1 + eps * k1) ** 2 + (r2 + eps * k2) ** 2 + (r3 + eps * k3) ** 2 / (1 - nu)) / G0n
    errs.append(np.max(np.abs(Gfull - 1 - eps * g(sN))))
ratio = errs[0] / errs[1]
report('G along the trapped front equals the toughness to first order (slope terms cancel)', errs[1] < 1e-5 and 50 < ratio < 200,
       'max |G/G0 - 1 - eps g| = %.2e at eps = 0.01, %.2e at eps = 0.001 (ratio %.0f; second order gives 100)' % (errs[0], errs[1], ratio))

# 7. maximum of kappa3 in closed form. kappa3/rho = a g + b S with S = 2 q sin(th)/(1 - 2 q cos(th) + q^2), th = 2 pi s;
#    d/dth = 0 gives b (1 + q^2) cos(th) - a (1 - q^2) sin(th) = 2 q b, i.e. R cos(th + phi) = 2 q b
a3, b3 = (2 + nu) / (2 - nu) * r3 / 2, 2 * (1 - nu) / (2 - nu) * r2
R, phi = math.hypot(b3 * (1 + q * q), a3 * (1 - q * q)), math.atan2(a3 * (1 - q * q), b3 * (1 + q * q))
def kap3(th):
    dn = 1 - 2 * q * math.cos(th) + q * q
    return rho * (a3 * ((1 - q * q) / dn - 1) + b3 * 2 * q * math.sin(th) / dn)
th_max = max((-phi + sgn * math.acos(2 * q * b3 / R) for sgn in (1.0, -1.0)), key=kap3)
s_cf = ((th_max / (2 * math.pi) + 0.5) % 1.0) - 0.5
report('maximum of kappa3: closed form against the reference search',
       abs(kap3(th_max) - ref['partB']['kappa3_max']) < 1e-10 and abs(s_cf - ref['partB']['s_kappa3_max']) < 1e-8,
       'kappa3_max %.12f at s = %.12f; reference differs by %.1e and %.1e' %
       (kap3(th_max), s_cf, abs(kap3(th_max) - ref['partB']['kappa3_max']), abs(s_cf - ref['partB']['s_kappa3_max'])))

# 8. loading level: G0 = Gamma0 at first order
report('K10_scaled', abs(1 / math.sqrt(G0n) - ref['partB']['K10_scaled']) < 1e-12, '%.12f' % (1 / math.sqrt(G0n)))
print('FAILED %d' % fails if fails else 'ALL CHECKS PASSED')
sys.exit(1 if fails else 0)
