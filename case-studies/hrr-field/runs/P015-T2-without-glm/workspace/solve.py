import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, brentq

np.seterr(all='ignore')

# ----------------------------------------------------------------------------
# HRR-type notch-tip eigenproblem (plane strain, pressure-insensitive power law)
#
# sigma_ij = C r^{-p} sig_i(θ);  Airy Φ = C r^{2-p} f(θ)
#   s_tt = (2-p)(1-p) f ;  s_rr = (2-p) f + f'' ;  s_rt = -(1-p) f'
#   D := s_rr - s_tt = f'' + p(2-p) f ;  T := s_rt ;  s_e^2 = (3/4)D^2 + 3T^2
# ε ∝ r^{-np} (plastic dominates) ; plane strain => s_zz=0 => σ_zz=(σ_rr+σ_θθ)/2
# Compatibility (u ∝ r^{1-np}; note factors: U = e/(1-np), V' = -(2-np)/(1-np) e):
#   (E^{n-1} D)'' + a (E^{n-1} D) - 4(1-np) (E^{n-1} T)' = 0 ,  a = np(2-np)
# BCs: flanks θ=±α traction-free => f(α)=f'(α)=0; symmetry => f'(0)=f'''(0)=0;
#      normalize f(0)=1 (amplitude free). Unknowns: p and q=f''(0).
# ----------------------------------------------------------------------------

def rhs_maker(p, n):
    c0 = p*(2.0-p); a = n*p*(2.0-n*p); m = n-1.0; c1 = 4.0*(1.0-n*p)
    def rhs(th, y):
        f, f1, f2, f3 = y
        D  = f2 + c0*f
        Dp = f3 + c0*f1
        T  = -(1.0-p)*f1; Tp = -(1.0-p)*f2; Tpp = -(1.0-p)*f3
        E2 = 0.75*D*D + 3.0*T*T
        if not np.isfinite(E2) or E2 < 1e-13:
            raise ZeroDivisionError('s_e -> 0 or overflow')
        if max(abs(f), abs(f1), abs(f2), abs(f3)) > 1e6:
            raise OverflowError('trajectory blows up (bad branch)')
        E  = np.sqrt(E2)
        Ep = (1.5*D*Dp + 6.0*T*Tp)/(2.0*E)
        num = m*((n-3.0)*Ep*Ep*D/E2 + D*(1.5*Dp*Dp + 6.0*Tp*Tp + 6.0*T*Tpp)/(2.0*E2)) \
              + 2.0*m*Ep*Dp/E + a*D - c1*(m*Ep*T/E + Tp)
        den = 1.0 + 0.75*m*D*D/E2
        Dpp = -num/den
        return [f1, f2, f3, Dpp - c0*f2]
    return rhs

def shoot(p, q, n, alpha, rtol=1e-11, atol=1e-13, dense=True):
    try:
        sol = solve_ivp(rhs_maker(p, n), (0.0, alpha), [1.0, 0.0, q, 0.0],
                        method='DOP853', rtol=rtol, atol=atol, dense_output=dense)
    except Exception:
        return None
    return sol if sol.success else None

def residual(p, q, n, alpha, rtol=1e-11, atol=1e-13):
    sol = shoot(p, q, n, alpha, rtol, atol, dense=False)
    if sol is None:
        return np.array([1e8, 1e8])
    return np.array([sol.y[0, -1], sol.y[1, -1]])

def find_roots(n, alpha, prange=(0.01, 0.249, 70), qrange=(-6, 2, 81), nseeds=30, rtol=1e-9, atol=1e-11):
    ps = np.linspace(*prange); qs = np.linspace(*qrange)
    R = np.array([[residual(p, q, n, alpha, rtol, atol) for q in qs] for p in ps])
    N = np.linalg.norm(R, axis=-1)
    roots = []
    for k in np.argsort(N.ravel())[:nseeds]:
        i, j = divmod(int(k), len(qs))
        if N[i, j] > 0.3:
            break
        try:
            s = least_squares(lambda x: residual(x[0], x[1], n, alpha, rtol, atol), [ps[i], qs[j]],
                              bounds=([1e-4, -40], [0.2499, 40]),
                              xtol=3e-16, ftol=3e-16, gtol=3e-16)
        except Exception:
            continue
        if s.cost < 1e-18 and np.all(np.isfinite(s.x)):
            p, q = float(s.x[0]), float(s.x[1])
            if 5e-4 < p < 0.2495 and not any(abs(p - r0) < 1e-6 for r0, _ in roots):
                roots.append((p, q))
    return sorted(roots, key=lambda t: -t[0])

# ----------------------------------------------------------------------------
# analytic elastic (n=1) Williams mode-I exponent for cross-check: p = 1 - λ,
# sin(2λα) + λ sin(2α) = 0
# ----------------------------------------------------------------------------
def williams_p(alpha):
    g = lambda lam: np.sin(2*lam*alpha) + lam*np.sin(2*alpha)
    lam = brentq(g, 1e-12, 1.0 - 1e-12)
    return 1.0 - lam

# ----------------------------------------------------------------------------
def report(n, alpha, p, q, label=''):
    sol = shoot(p, q, n, alpha, rtol=1e-12, atol=1e-14)
    th = np.linspace(0, alpha, 4001)
    Y = sol.sol(th)
    f, f1, f2 = Y[0], Y[1], Y[2]
    c0 = p*(2-p); a = n*p*(2-n*p); c1 = 4*(1-n*p)
    srr = (2-p)*f + f2
    stt = (2-p)*(1-p)*f
    srt = -(1-p)*f1
    D = f2 + c0*f
    E = np.sqrt(0.75*D*D + 3*srt*srt)
    G = E**(n-1)
    GD = G*D; GT = G*srt
    L = np.gradient(np.gradient(GD, th), th) + a*GD - c1*np.gradient(GT, th)
    scale = max(np.max(np.abs(a*GD)), np.max(np.abs(np.gradient(GT, th))))
    eq1 = np.max(np.abs(np.gradient(srt, th) - (stt - (1-p)*srr)))
    eq2 = np.max(np.abs(np.gradient(stt, th) - (p-2)*srt))
    R0 = ((2-p) + q)/((2-p)*(1-p))
    out = dict(
        label=label, p=p, q=q, R=R0,
        bc=np.array([sol.y[0, -1], sol.y[1, -1]]),
        minE2=float(np.min(E**2)),
        ode_rel=float(np.max(np.abs(L))/scale),
        eq=float(max(eq1, eq2)),
        stt0=float(stt[0]), srr0=float(srr[0]), szz0=float(0.5*(srr[0]+stt[0])),
        se_theta={float(np.round(t, 4)): float(np.interp(t, th, E)) for t in [0, alpha/4, alpha/2, 3*alpha/4]},
    )
    print(f"[{label}] p={p:.8f} q={q:.6f} R={R0:.6f} BC_res={out['bc']} minE2={out['minE2']:.2e} "
          f"ODE_rel={out['ode_rel']:.2e} eq_res={out['eq']:.2e}")
    return out

# ----------------------------------------------------------------------------
if __name__ == '__main__':
    tw = 2*np.pi/3  # half material angle α (flanks at ±2π/3, material angle 240°)

    print('--- validations: crack (material angle 360deg) must give p = 1/(n+1) ---')
    for n, pt in [(1.0, 0.5), (3.0, 0.25), (4.0, 0.2), (9.0, 0.1)]:
        r = find_roots(n, np.pi, prange=(max(0.3*pt, pt-0.15), min(pt+0.15, 0.6), 31),
                       qrange=(-4, 1, 41), nseeds=12)
        best = [(round(p, 8), round(q, 5), round(abs(p-pt), 10)) for p, q in r]
        print(f"n={n:.0f} crack roots (target p={pt}):", best)

    print('--- main: n=4, material angle 240 deg ---')
    roots = find_roots(4.0, tw)
    print("all roots (p,q):", [(round(p, 8), round(q, 5)) for p, q in roots])
    p, q = roots[0]
    main = report(4.0, tw, p, q, 'main')
