"""Local grading audit: use only delivered definitions, without running solve.py."""
from pathlib import Path
import json

namespace = {}
source = Path('scratch/reproduction/solve.py').read_text()
exec(compile(source.split('# --- convergence')[0], 'submitted_definitions', 'exec'), namespace)
mp = namespace['mp']
sp = namespace['sp']
x = namespace['x']

def exact_basis(N):
    result = [[mp.mpf(0), mp.mpf(3), mp.mpf(-3), mp.mpf(1)]]
    for j in range(N + 1):
        polynomial = sp.Poly(sp.expand(x * (1-x)**3 * sp.legendre(j, 2*x-1)), x)
        result.append([mp.mpf(int(polynomial.nth(k))) for k in range(polynomial.degree()+1)])
    return result

def endpoint_values(c):
    return [namespace['peval'](namespace['pderiv'](c, k), mp.mpf(1)) for k in (0, 1, 2)]

print('Submitted expanded Legendre basis: endpoint defects (should all vanish except phi0(1)=1)')
for N in (8, 12, 16, 20, 24, 28):
    approximate = namespace['make_basis'](N)
    exact = exact_basis(N)
    errors = [max(abs(mp.mpf(float(v))-w) for v, w in zip(c, d)) for c, d in zip(approximate, exact)]
    defect = [max(abs(endpoint_values(c)[k]) for c in approximate[1:]) for k in range(3)]
    print('N=%d max_coefficient_error=%s max_endpoint_defects=%s' %
          (N, mp.nstr(max(errors), 10), [mp.nstr(v, 10) for v in defect]), flush=True)

# The purported alternative basis has a coefficient-order bug before any collocation.
import numpy as np
w = np.array([1.0, -1.0])
one_plus = np.polyadd(np.array([1.0]), w)
one_plus = np.polyadd(one_plus, np.convolve(w, w))
wrong_phi0 = np.convolve(np.array([0.0, 1.0]), one_plus)
print('explore26 phi0 ascending coefficients:', wrong_phi0.tolist())
print('explore26 phi0 at (1), first and second derivatives:', [mp.nstr(v, 10) for v in endpoint_values(wrong_phi0)])

refs = {'lambda1': 11.74265, 'lambda2': 264.6348, 'mode1_shape_half': 0.8020}
claims = json.loads(Path('submission/results.json').read_text())
print('Relative numerical errors:', {k: (claims[k]/v-1)*100 for k, v in refs.items()})

# The submitted Q=80 rule is not exact for the rational reflected-kernel term.
for xv in ('0.4', '0.7', '0.004657026981834667'):
    xx = mp.mpf(xv)
    c = [mp.mpf(0), mp.mpf(3), mp.mpf(-3), mp.mpf(1)]
    submitted = namespace['Nphi'](xx, c)
    # Exact integral of S for a constant third derivative (phi0'''=6).
    integral_S = mp.log((1+xx)/(1-xx)) + 6*xx/(1+xx) - 2*xx**2/(1+xx)**2 - 4
    actual = 6*integral_S
    print('phi0 operator x=%s submitted=%s analytic=%s error=%s' %
          (xv, mp.nstr(submitted, 18), mp.nstr(actual, 18), mp.nstr(submitted-actual, 12)), flush=True)
