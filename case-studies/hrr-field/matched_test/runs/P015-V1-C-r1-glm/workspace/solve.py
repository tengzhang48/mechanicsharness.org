import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

TAU0, GAM0, JGIVEN = 80.0, 0.002, 0.120
N = 3.5


def rhs_factory(n):
    a = 1.0 / n - 1.0

    def rhs(theta, y, s):
        f, fp = y
        G2 = s * s * f * f + fp * fp
        c0 = s * (1.0 + (s - 1.0) / n)
        fpp = -f * (a * s * s * fp * fp + c0 * G2) / (s * s * f * f + fp * fp / n)
        return (fp, fpp)

    return rhs


def ivp(n, s, dense):
    return solve_ivp(rhs_factory(n), (0.0, np.pi), (0.0, 1.0), args=(s,),
                     method='DOP853', rtol=1e-12, atol=1e-14, dense_output=dense)


def find_roots(n, lo=0.03, hi=0.60, m=58):
    def resid(s):
        return ivp(n, s, False).y[1, -1]
    grid = np.linspace(lo, hi, m)
    vals = np.array([resid(s) for s in grid])
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(resid, grid[i], grid[i + 1], xtol=1e-15, rtol=1e-15))
    return roots


def angular(n, s, m=200001):
    sol = ivp(n, s, True)
    th = np.linspace(0.0, np.pi, m)
    f, fp = sol.sol(th)
    G = np.sqrt(s * s * f * f + fp * fp)
    p = (n + 1.0) / n
    integrand = (n / (n + 1.0)) * G**p * np.cos(th) \
        - s * f * G**(1.0 / n - 1.0) * (s * f * np.cos(th) - fp * np.sin(th))
    Iang = 2.0 * np.trapezoid(integrand, th)
    return sol, Iang, f[-1], fp[-1]


def amplitude(n, Iang, J):
    return (J / (TAU0 * GAM0 ** (-1.0 / n) * Iang)) ** (n / (n + 1.0))


def field(x, y, n, s, sol, A):
    x = np.atleast_1d(np.asarray(x, float))
    y = np.atleast_1d(np.asarray(y, float))
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    sgn = np.where(th >= 0.0, 1.0, -1.0)
    f, fp = sol.sol(np.abs(th))
    f = sgn * f
    rp = r ** (s - 1.0)
    gr = A * s * rp * f          # gamma . e_r
    gth = A * rp * fp            # gamma . e_theta
    gx = gr * np.cos(th) - gth * np.sin(th)
    gy = gr * np.sin(th) + gth * np.cos(th)
    g = np.hypot(gx, gy)
    tau = TAU0 * (g / GAM0) ** (1.0 / n)
    tx = tau * gx / g
    ty = tau * gy / g
    W = (n / (n + 1.0)) * TAU0 * GAM0 ** (-1.0 / n) * g ** ((n + 1.0) / n)
    return dict(gx=gx, gy=gy, g=g, gr=gr, gth=gth, tx=tx, ty=ty, tau=tau, W=W)


def J_rect(n, s, sol, A, a=0.10, m=40001):
    total = 0.0
    dx = 2 * a / (m - 1)

    def add(xs, ys, Nx, Ny):
        F = field(xs, ys, n, s, sol, A)
        tz = F['tx'] * Nx + F['ty'] * Ny
        return np.trapezoid(F['W'] * Nx - tz * F['gx'], dx=dx)

    t = np.linspace(-a, a, m)
    total += add(t, np.full(m, -a), 0.0, -1.0)   # bottom, outward N=(0,-1)
    total += add(np.full(m, a), t, 1.0, 0.0)     # right
    total += add(t, np.full(m, a), 0.0, 1.0)     # top
    total += add(np.full(m, -a), t, -1.0, 0.0)   # left
    return total


def J_arc(n, s, sol, A, r, m=400001):
    th = np.linspace(-np.pi, np.pi, m)
    F = field(r * np.cos(th), r * np.sin(th), n, s, sol, A)
    tz = F['tx'] * np.cos(th) + F['ty'] * np.sin(th)   # tau . e_r
    return np.trapezoid(F['W'] * np.cos(th) - tz * F['gx'], th) * r


def main():
    print('=== n = 3.5 main case ===')
    roots = find_roots(N)
    s_th = 1.0 / (N + 1.0)
    print('theory exponent 1/(n+1) =', s_th)
    print('eigenvalue roots found  =', roots)
    s = min(roots, key=lambda r: abs(r - s_th))
    sol, Iang, fpi, fppi = angular(N, s)
    A = amplitude(N, Iang, JGIVEN)
    print(f's* = {s:.15f}  (theory {s_th:.15f}, diff {abs(s-s_th):.2e})')
    print(f'f(pi) = {fpi:.12f}   fprime(pi) = {fppi:.3e}   I_ang = {Iang:.12f}')
    print(f'A = {A:.12e} mm^(1-s)')

    r1, r2 = 0.0200, 0.0800
    dw1 = 2 * A * r1**s * fpi
    dw2 = 2 * A * r2**s * fpi
    print(f'\n1) face slip Dw(r): {dw1:.6e} mm at r=0.02 ; {dw2:.6e} mm at r=0.08')
    print(f'   ratio Dw(0.08)/Dw(0.02) = {dw2/dw1:.6f} (theory (4)^s = {4**s:.6f})')

    th60 = np.pi / 3
    F = field(r1 * np.cos(th60), r1 * np.sin(th60), N, s, sol, A)
    trz = F['tx'][0] * np.cos(th60) + F['ty'][0] * np.sin(th60)
    tthz = -F['tx'][0] * np.sin(th60) + F['ty'][0] * np.cos(th60)
    print(f'\n2) at r=0.02, theta=60deg: tau_rz = {trz:.6f} MPa ; tau_thz = {tthz:.6f} MPa')
    print(f'   (gamma there = {F["g"][0]:.6e}, tau = {F["tau"][0]:.6f} MPa)')

    gtarget = 0.0100
    r0 = (A * 1.0 / gtarget) ** (1.0 / (1.0 - s))          # G(0)=1
    G90 = float(np.sqrt((s * sol.sol(np.pi / 2)[0])**2 + (sol.sol(np.pi / 2)[1])**2))
    Gpi = s * fpi                                          # f'(pi)=0
    r90 = (A * G90 / gtarget) ** (1.0 / (1.0 - s))
    rpi = (A * Gpi / gtarget) ** (1.0 / (1.0 - s))
    print(f'\n3) |gamma|=0.01 radii: theta=0: {r0:.6f} ; 90deg: {r90:.6f} ; 180deg(face): {rpi:.6f} mm')
    print(f'   G(0)=1, G(90)={G90:.9f}, G(180)={Gpi:.9f}')

    print('\n=== verification: J path checks (input J = 0.120) ===')
    print(f'J on circle r=0.02 : {J_arc(N, s, sol, A, 0.02):.9f}')
    print(f'J on circle r=0.15 : {J_arc(N, s, sol, A, 0.15):.9f}')
    print(f'J on square a=0.10 : {J_rect(N, s, sol, A, 0.10):.9f}')
    print(f'J on square a=0.25 : {J_rect(N, s, sol, A, 0.25):.9f}')

    print('\n=== verification: linear limit n=1 vs closed-form K-field ===')
    roots1 = find_roots(1.0)
    print('n=1 roots:', roots1)
    s1 = min(roots1, key=lambda r: abs(r - 0.5))
    sol1, Iang1, fpi1, fppi1 = angular(1.0, s1)
    print(f's1 = {s1:.15f} (theory 0.5) ; I_ang = {Iang1:.12f} (theory pi = {np.pi:.12f}) ; f(pi) = {fpi1:.12f} (theory 2)')
    mu = TAU0 / GAM0
    K = 40.0
    J1 = K**2 / (2 * mu)
    A1 = amplitude(1.0, Iang1, J1)
    dw_num = 2 * A1 * np.sqrt(r1) * fpi1
    dw_ex = K * np.sqrt(8 * r1 / np.pi) / mu
    print(f'with K=40 MPa sqrt(mm): J = K^2/(2 mu) = {J1:.9f} N/mm')
    print(f'n=1 pipeline Dw(0.02) = {dw_num:.9e} ; closed form K sqrt(8r/pi)/mu = {dw_ex:.9e}')
    print(f'n=1 J_rect check: {J_rect(1.0, s1, sol1, A1, 0.1):.9f} vs {J1:.9f}')

    print('\n=== face residual ===')
    th = np.pi - 1e-9
    f, fp = sol.sol(th)
    print(f'|sigma_thz/sigma_rz| on face ~ |fprime(pi)|/(s f(pi)) = {abs(fppi)/(s*abs(fpi)):.3e}')


main()
