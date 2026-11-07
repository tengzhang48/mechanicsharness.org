"""Audit the submitted spectrum, then use exact coefficients and integrated kernel."""
from pathlib import Path
import mpmath as mp
import sympy as sp

ns = {}
source = Path('scratch/reproduction/solve.py').read_text()
exec(compile(source.split('# --- convergence')[0], 'submitted_definitions', 'exec'), ns)
mp.mp.dps = 70
x = ns['x']

def exact_basis(N):
    b = [[mp.mpf(0), mp.mpf(3), mp.mpf(-3), mp.mpf(1)]]
    for j in range(N+1):
        p = sp.Poly(sp.expand(x*(1-x)**3*sp.legendre(j, 2*x-1)), x)
        b.append([mp.mpf(int(p.nth(k))) for k in range(p.degree()+1)])
    return b

def monomial_operator_values(xv, degree):
    # A_n,B_n,C_n are integrals of xi^n/(x+xi)^1,^2,^3.
    a = mp.log((1+xv)/xv)
    b = 1/xv - 1/(1+xv)
    c = (1/xv**2 - 1/(1+xv)**2)/2
    H = mp.log(xv/(1-xv))
    result = [mp.mpf(0)]*3
    for n in range(degree-2):
        if n:
            a, b, c = 1/mp.mpf(n)-xv*a, a-xv*b, b-xv*c
        hilbert = xv**n*H - sum(xv**k/mp.mpf(n-k) for k in range(n))
        val = hilbert+a-6*xv*b+4*xv**2*c
        p = n+3
        result.append(p*(p-1)*(p-2)*val)
    return result

def spectrum(N, exact=False):
    basis = exact_basis(N) if exact else ns['make_basis'](N)
    m = len(basis)
    # Same submitted nodes for direct comparison; exactly generated for exact runs.
    if exact:
        nodes = [mp.mpf('0.5')*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
    else:
        nodes = [mp.mpf(str(0.5*(1-ns['np'].cos(ns['np'].pi*(i+1)/(m+1))))) for i in range(m)]
    O = mp.matrix(m); V = mp.matrix(m)
    for i, xi in enumerate(nodes):
        if exact:
            mon = monomial_operator_values(xi, N+4)
        for k, p in enumerate(basis):
            V[i,k] = ns['peval'](p, xi)
            O[i,k] = sum(p[j]*mon[j] for j in range(len(p))) if exact else ns['Nphi'](xi,p)
    L = mp.inverse(V)*O
    ev = mp.eig(L, left=False, right=False)
    nonreal = [e for e in ev if abs(mp.im(e)) > mp.mpf('1e-12')*max(1,abs(e))]
    nondecaying = [e for e in ev if mp.re(e) >= 0]
    pos = sorted([-mp.re(e) for e in ev if abs(mp.im(e)) < mp.mpf('1e-12')*max(1,abs(e)) and mp.re(e)<0])
    return pos, nonreal, nondecaying

for N, exact in ((20,False),(20,True),(24,True),(28,True),(40,True)):
    rates, nonreal, nondecaying = spectrum(N,exact)
    print('N=%d exact_coefficients_and_integrals=%s first_rates=%s nonreal=%d nondecaying=%d' %
          (N, exact, [mp.nstr(v,14) for v in rates[:3]], len(nonreal),len(nondecaying)),flush=True)
    if not exact:
        print('Submitted discarded eigenvalues:', [mp.nstr(e,15) for e in nonreal+nondecaying],flush=True)

# Check the derived integrated monomial kernel against direct PV subtraction.
for xv in ('0.4','0.7'):
    xx = mp.mpf(xv)
    val = monomial_operator_values(xx,8)
    for p in (3,4,5,8):
        q = ns['Nphi'](xx,[mp.mpf(0)]*p+[mp.mpf(1)])
        print('monomial cross-check x=%s p=%d error=%s' % (xv,p,mp.nstr(q-val[p],6)),flush=True)
