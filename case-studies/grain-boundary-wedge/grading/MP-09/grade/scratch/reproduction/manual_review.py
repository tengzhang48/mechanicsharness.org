"""Targeted grader diagnostics; all inputs are scratch copies of the submission."""
import contextlib
import io
import json
import numpy as np
from scipy.linalg import eig
import solve
import sblib as sb
from gblib import T_moment

original = json.load(open('../../submission/results.json'))
rerun = json.load(open('results.json'))
print('Reproduced JSON identical:', original == rerun)
print('Required value errors relative to key:')
for field, reference in [('lambda1', 11.74265), ('lambda2', 264.6348),
                         ('mode1_shape_half', .802), ('K_inf_norm', 1.1215)]:
    print(field, 'value=', rerun[field], 'relative_error=', rerun[field]/reference-1)

print('\nFD manufactured-function errors, shrinking versus fixed interior:')
for option, builder in [('A', solve.build_FD_cellcentered),
                        ('B', solve.build_FD_nodebased)]:
    for n in [64, 256, 1024]:
        L, x, _ = builder(n)
        f = x*x*(1-x)**3
        exact = np.array([-18*T_moment(t,0)+72*T_moment(t,1)-60*T_moment(t,2) for t in x])
        err = np.abs(L@f-exact)
        j = 4+np.argmax(err[4:-4])
        fixed = (x>=.05)&(x<=.95)
        print(option, n, 'four-cell-trimmed_max=', err[j], 'location=', x[j],
              'fixed_[.05,.95]_max=', np.max(err[fixed]))

print('\nSpectral tolerance/resolution sensitivity:')
cache = {}
for n in [60, 80, 110]:
    ss, Th, F = sb.build_s_matrices(n)
    U, singular, Vt = np.linalg.svd(F)
    cache[n] = (ss, Th, F, U, singular, Vt)
    for tol in [1e-11, 1e-12, 1e-13]:
        rank = np.count_nonzero(singular > tol*singular[0])
        Ur, Vr = U[:,:rank], Vt[:rank].T
        w, vectors = eig(Ur.T@Th@Vr, Ur.T@F@Vr)
        ok = np.isfinite(w)&(np.abs(w.imag)<1e-7*np.maximum(1,np.abs(w.real)))&(w.real<0)
        rates = np.sort(-w[ok].real)
        print('N=',n,'tol=',tol,'rank=',rank,'lambda1=',rates[0],'lambda2=',rates[1])

lambda1, lambda2, half, a = sb.mode1_spectral_s(110)
ks = np.arange(5,5+len(a))
print('\nTip expansion coefficient a5=',a[0],'; a7=',a[2])
for r in [1e-3, 1e-4, 1e-5, 1e-6]:
    s=np.sqrt(r)
    delta=np.sum(a*s**ks)
    third=np.sum(a*sb.c_coeff(ks)*s**(ks-6))
    print('r=',r,'(f(1-r)-1)/r^2.5=',delta/r**2.5,'sqrt(r)*fthird=',np.sqrt(r)*third)

print('\nTrue non-normal spectral projection of uniform initial data:')
# Use a small independent FD grid to verify that initial/boundary point
# incompatibility does not prevent a well-defined biorthogonal projection.
L, x, _=solve.build_FD_cellcentered(128)
w, left, right=eig(L,left=True,right=True)
i=np.argmax(w.real)
v=right[:,i].real
v/=v[-1]
ell=left[:,i].real
coefficient=(ell@np.ones(len(x)))/(ell@v)
print('FD128 lambda1=',-w[i].real,'biorthogonal c1=',coefficient)
print('Tip BC counterexample: s^4 has fs(0)=fss(0)=fsss(0)=0, but fxx(1)=2.')

print('\nCrack matrix conditioning:')
for m in [10,14,18,24]:
    points=sb.cheb_points_s(m)
    A=np.array([2*sb.R_all(s,2*m-2)[2*np.arange(m)] for s in points])
    c=np.linalg.solve(A,np.ones(m))
    print('M=',m,'cond(A)=',np.linalg.cond(A),'maxabs(c)=',np.max(np.abs(c)))

print('\nNear-tip COD-K coefficient check:')
crack=sb.crack_K_and_COD(24)
c=crack['c']
expected=8*crack['K_tip']/np.sqrt(2)
for r in [1e-2,1e-4,1e-6,1e-8]:
    # COD(r) = -8*pi sum c_m*r^(m+1/2)/(2*m+1).
    coefficient=-8*np.pi*np.sum(c*r**np.arange(len(c))/(2*np.arange(len(c))+1))
    print('r=',r,'COD/sqrt(r)=',coefficient,'8*Knorm/sqrt(2)=',expected)
print('Identity at r=0: -8*pi*c0=',-8*np.pi*c[0],'; same c0 determines K_tip.')

print('\nOutside stress comparison with an intentionally wrong crack density:')
wrong=np.array([-0.1])
outside=[]
for e in [1e-4,1e-6,1e-8]:
    outside.append(np.sqrt(2*e)*(1-sb.T_outside(wrong,1+e)))
print('Wrong c0=',wrong[0],'K_tip=',-np.pi*np.sqrt(2)*wrong[0],
      'K_out_samples=',outside,'yet boundary traction equation is not satisfied.')

print('\nStored Richardson interval versus spectral value:')
print('11.74267930 <= lambda1 <= 11.74268899:',11.74267930<=lambda1<=11.74268899)
