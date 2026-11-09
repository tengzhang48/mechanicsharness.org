import numpy as np
from scipy.optimize import least_squares
import solve as S

p, q, n, alpha = 0.15617526, -1.263079, 4.0, 2*np.pi/3

# 1) ODE residual via FD vs grid density; also profile of trajectory
for Npts in [401, 2001, 8001, 32001]:
    sol = S.shoot(p, q, n, alpha, rtol=1e-13, atol=1e-15)
    th = np.linspace(0, alpha, Npts)
    Y = sol.sol(th)
    c0 = p*(2-p); a = n*p*(2-n*p); c1 = 4*(1-n*p)
    D = Y[2] + c0*Y[0]
    T = -(1-p)*Y[1]
    E = np.sqrt(0.75*D*D + 3*T*T)
    G = E**(n-1); GD = G*D; GT = G*T
    L = np.gradient(np.gradient(GD, th), th) + a*GD - c1*np.gradient(GT, th)
    scale = max(np.max(np.abs(a*GD)), np.max(np.abs(c1*np.gradient(GT, th))),
                np.max(np.abs(np.gradient(np.gradient(GD, th), th))))
    print(f"Npts={Npts:6d}  ODE_rel={np.max(np.abs(L))/scale:.3e}  "
          f"max|f|={np.max(np.abs(Y[0])):.3f} max|f3|={np.max(np.abs(Y[3])):.3f} nsteps={sol.t.size}")

# 2) polish the root at tight ODE tolerance
solref = lambda x: S.residual(x[0], x[1], n, alpha, rtol=1e-13, atol=1e-15)
s = least_squares(solref, [p, q], bounds=([1e-4, -10], [0.249, 10]),
                  xtol=1e-15, ftol=1e-15, gtol=1e-15)
p2, q2 = s.x
R2 = ((2-p2)+q2)/((2-p2)*(1-p2))
print(f"polished: p={p2:.10f} q={q2:.8f} R={R2:.8f} res={s.fun} cost={s.cost:.2e}")
# sensitivity: perturb p by +-3e-4 (tolerance scale)
for dp in [-3e-4, 3e-4]:
    sp = least_squares(solref, [p2+dp, q2], bounds=([1e-4, -10], [0.249, 10]))
    print(f"  seed p+{dp:+.0e}: -> p={sp.x[0]:.10f} (same root? {abs(sp.x[0]-p2)<1e-8})")

# 3) missing crack validations with relaxed cost threshold
def find_roots2(n, alpha, prange, qrange, nseeds=12):
    ps = np.linspace(*prange); qs = np.linspace(*qrange)
    R = np.array([[S.residual(pp, qq, n, alpha, 1e-10, 1e-12) for qq in qs] for pp in ps])
    N = np.linalg.norm(R, axis=-1)
    roots = []
    for k in np.argsort(N.ravel())[:nseeds]:
        i, j = divmod(int(k), len(qs))
        if N[i, j] > 0.5: break
        try:
            s = least_squares(lambda x: S.residual(x[0], x[1], n, alpha, 1e-12, 1e-14),
                              [ps[i], qs[j]], bounds=([1e-4, -40], [0.9, 40]),
                              xtol=1e-15, ftol=1e-15, gtol=1e-15)
        except Exception:
            continue
        if s.cost < 1e-14 and abs(s.x[0]-prange[0]) > 1e-3:
            if not any(abs(s.x[0]-r0) < 1e-6 for r0, _ in roots):
                roots.append((float(s.x[0]), float(s.x[1])))
    return roots
for n, pt in [(1.0, 0.5), (3.0, 0.25)]:
    r = find_roots2(n, np.pi, (pt-0.12, pt+0.12, 25), (-3, 1, 41))
    print(f"n={n:.0f} crack roots (target {pt}):", [(round(a, 9), round(b, 5)) for a, b in r])

# 4) elastic-plausibility anchor: elastic Williams mode-I 240deg wedge R at theta=0
alpha_w = 2*np.pi/3
from scipy.optimize import brentq
gl = lambda lam: np.sin(2*lam*alpha_w) + lam*np.sin(2*alpha_w)
lam = brentq(gl, 1e-9, 1-1e-9); pw = 1-lam
c0e = pw*(2-pw); ae = pw*(2-pw)  # n=1: a=np(2-np)=p(2-p), c1=4(1-p)
# k^4 - [c0+a+4(1-p)]k^2 + a*c0 = 0
b2 = c0e + ae + 4*(1-pw)
ksq = np.roots([1, 0, -b2, 0, ae*c0e])
ks = np.sqrt([k for k in ksq if abs(k.imag) < 1e-9 and k.real > 0])
k1, k2 = np.sort(ks.real)
# eigenfunction: A cos k1 th + B cos k2 th with BCs at alpha
M = np.array([[np.cos(k1*alpha_w), np.cos(k2*alpha_w)],
              [k1*np.sin(k1*alpha_w), k2*np.sin(k2*alpha_w)]])
# null vector of M (A,B): f(alpha)=0, f'(alpha)=0
_, _, Vt = np.linalg.svd(M)
A, B = Vt[-1]
# normalize f(0)=1
A, B = A/(A+B), B/(A+B)
Re = ((2-pw)*(A+B) - (A*k1**2 + B*k2**2))/((2-pw)*(1-pw)*(A+B))
print(f"elastic 240deg wedge: lambda={lam:.6f} p={pw:.6f} k1={k1:.4f} k2={k2:.4f} R_elastic={Re:.4f}")
