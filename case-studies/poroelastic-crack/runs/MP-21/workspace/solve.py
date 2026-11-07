#!/usr/bin/env python3
"""
MP-21: a suddenly pressurized crack in a fluid-saturated porous solid.

Plane-strain Biot poroelasticity (Rice-Cleary form), Griffith crack -L<=x<=L.
Loading 1: sigma_n = -1, p = 0 on the faces.  Loading 2: sigma_n = 0, p = 1 on faces.

Analytic solution used (derived in ANSWER.md):
  Loading 1:
    t=0+  undrained elastic response (nu -> nu_u): K = sqrt(pi L),
          V = pi (1-nu_u) L^2/G;  induced p on faces = 2(nu_u-nu)/(alpha(1-2nu)).
    t->inf drained elastic response (p = 0 everywhere): K = sqrt(pi L),
          V = pi (1-nu) L^2/G.
    small times: drainage boundary layer on the faces (local oedometer problem,
    effective diffusivity exactly c) gives
          V ~ V1_0 + [8(nu_u-nu)/(sqrt(pi)(1-nu))] * sqrt(ct/L^2) * L^2/G.
  Loading 2:
    t=0+ nothing has diffused: K2 = 0.
    t->inf: 2D recurrence -> p -> 1 pointwise everywhere (no decaying steady state exists;
    the harmonic function 1-Re[z/sqrt(z^2-L^2)] with p->0 at infinity is unbounded at the
    tips and is NOT the parabolic limit - checked numerically). With p->1 and total stress
    -> 0 the final mechanical state is uniform swelling e_ij = alpha(1-2nu)/(2G) delta_ij,
    u = beta_sw*(x,y), which has zero relative displacement across the crack:
          K2_inf = 0, V2_inf = 0.
    (Transiently, V2(t) < 0 and K2(t) < 0: infiltration swells the matrix near the faces
    and pushes them toward each other.)

Checks:
 (1) elastic crack: independent Chebyshev-modal distributed-dislocation solve and
     Westergaard closed form vs analytic K, V, face stresses (drained nu and undrained nu_u).
 (2) local 1D poroelastic boundary layer: implicit FD vs erf profile and vs the
     surface-displacement formula that produces V1_short_coeff.
 (3) loading-2 long-time state: (a) NEGATIVE check: FD Laplace box solution vs the
     unbounded harmonic field 1-Re[z/sqrt(z^2-L^2)] (they differ -> that field rejected);
     (b) parabolic FD showing p -> 1 pointwise.
Writes results.json; prints every reported value.
"""
import json
import numpy as np
from scipy.special import erf
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

if not hasattr(np, "trapezoid"):                     # numpy < 2.0
    np.trapezoid = np.trapz

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, L, P0 = 1.0, 1.0, 1.0
eta = (nu_u - nu) / (alpha * (1 - 2 * nu))          # Skempton: p = -eta*sigma_kk (undrained)
beta_sw = alpha * (1 - 2 * nu) / (2 * G)            # e_ij = beta_sw * p * delta_ij at sigma=0

print("=" * 72)
print("MP-21  suddenly pressurized crack in a poroelastic solid")
print(f"nu={nu}  nu_u={nu_u}  alpha={alpha}  G={G}  L={L}  (per unit loading magnitude)")
print(f"eta=(nu_u-nu)/(alpha(1-2nu)) = {eta:.9f}   beta_sw=alpha(1-2nu)/(2G) = {beta_sw:.9f}")
print("=" * 72)

# --------------------------------------------------------------------------
# analytic results (general formulas evaluated for the parameters above)
# --------------------------------------------------------------------------
K1_0 = P0 * np.sqrt(np.pi * L)                      # per unit load
K1_inf = P0 * np.sqrt(np.pi * L)
V1_0 = np.pi * (1 - nu_u) * P0 * L**2 / G
V1_inf = np.pi * (1 - nu) * P0 * L**2 / G
V1_short_coeff = 8 * (nu_u - nu) / (np.sqrt(np.pi) * (1 - nu))
p_faces_0 = 2 * (nu_u - nu) / (alpha * (1 - 2 * nu))
K2_0 = 0.0
K2_inf = 0.0
V2_inf = 0.0


# ==========================================================================
# CHECK 1: elastic Griffith crack, independent numerical solutions
# ==========================================================================
def check1():
    print("\n--- CHECK 1: elastic crack (DD collocation vs Westergaard vs analytic) ---")

    # ---------- (a) Westergaard closed form ----------
    def west(z):
        s = np.sqrt(z * z - L * L + 0j)
        Z = P0 * (z / s - 1.0)
        dZ = -P0 * L**2 / (z * z - L * L + 0j) ** 1.5
        return Z, dZ

    def sig(x, y):
        Z, dZ = west(x + 1j * y)
        sxx = Z.real - y * dZ.imag
        syy = Z.real + y * dZ.imag
        txy = -y * dZ.real
        return sxx, syy, txy

    # equilibrium residual by finite differences at an interior point
    x0, y0, d = 0.37 * L, 0.53 * L, 1e-5
    dxsxx = (sig(x0 + d, y0)[0] - sig(x0 - d, y0)[0]) / (2 * d)
    dytxy = (sig(x0, y0 + d)[2] - sig(x0, y0 - d)[2]) / (2 * d)
    print(f"equilibrium d(sigma_xx)/dx + d(tau_xy)/dy at (0.37L,0.53L) = {dxsxx + dytxy:.3e}"
          f"   (sigma_xx={sig(x0,y0)[0]:.6f}, sigma_yy={sig(x0,y0)[1]:.6f})")

    # stresses on the crack faces (slightly above y=0)
    xf = np.array([-0.9, -0.5, 0.0, 0.5, 0.9]) * L
    sxxf = sig(xf, 1e-4)[0]
    syyf = sig(xf, 1e-4)[1]
    print("faces |x|/L =", np.round(xf / L, 2))
    print("  sigma_xx(face) =", np.round(sxxf, 6), " (expect -1)")
    print("  sigma_yy(face) =", np.round(syyf, 6), " (expect -1)")
    print(f"  -> sigma_kk(face) = -2  =>  p_u(face) = 2*eta = {2*eta:.6f}  (p_faces_0)")

    # K from the near-tip stress ahead of the tip (approaches from below by
    # P0*sqrt(2 pi r), the O(sqrt(r)) term of the exact field)
    for r in [1e-2, 1e-3, 1e-4]:
        k_raw = sig(L + r, 0.0)[1] * np.sqrt(2 * np.pi * r)
        print(f"  K from sigma_yy(L+{r:.0e},0): raw = {k_raw:.6f}, bias-corrected "
              f"= {k_raw + P0*np.sqrt(2*np.pi*r):.8f}  (expect {P0*np.sqrt(np.pi*L):.8f})")

    # crack face displacement from -int e_yy dy (drained), substitution y=s^2,
    # plus analytic far tail: e_yy -> (nu-3/2)*P0*L^2/(2G y^2) on the column x=0
    s = np.linspace(0, 8, 200001)
    yv = s * s
    for xq in [0.0, 0.5, 0.9]:
        sxx, syy, _ = sig(xq, yv)
        eyy = (syy - nu * (sxx + syy)) / (2 * G)
        tail = (nu - 1.5) * P0 * L**2 / (2 * G * yv[-1])
        uy = -np.trapezoid(2 * s * eyy, s) - tail  # dy = 2 s ds
        print(f"  u_y(x={xq}L) = {uy:.8f}  expect (1-nu)*sqrt(L^2-x^2)/G = "
              f"{(1-nu)*np.sqrt(L**2-xq**2):.8f}")

    # ---------- (b) distributed dislocation collocation solve ----------
    # unknown phi(xi_j) at Gauss-Chebyshev nodes, b = phi/sqrt(L^2-xi^2);
    # kernel: sigma_yy(x,0) = G/(2 pi (1-nu)) PV int b(xi)/(x-xi) dxi,  b = phi/sqrt(L^2-xi^2)
    # Using PV int w(t) T_n(t)/(t-a) dt = pi U_{n-1}(a)  (|a|<1), the SIE becomes
    # sum_n a_n U_{n-1}(x/L) = 2(1-nu) P(x) L/G  with phi = sum a_n T_n (odd n, symmetric load)
    def chebT_rec(nmax, x):
        T = np.zeros((nmax + 1, x.size))
        T[0] = 1.0
        if nmax >= 1:
            T[1] = x
        for k in range(2, nmax + 1):
            T[k] = 2 * x * T[k - 1] - T[k - 2]
        return T

    def chebU_rec(nmax, x):
        U = np.zeros((nmax + 1, x.size))
        U[0] = 1.0
        if nmax >= 1:
            U[1] = 2 * x
        for k in range(2, nmax + 1):
            U[k] = 2 * x * U[k - 1] - U[k - 2]
        return U

    # self-test of the PV identity by dense integration with singularity subtraction
    a_test = 0.3
    t = np.linspace(-1 + 1e-6, 1 - 1e-6, 400001)
    w = 1 / np.sqrt(1 - t * t)
    for n_test in [1, 2, 4]:
        Tn = chebT_rec(n_test, t)[n_test]
        Tn_a = chebT_rec(n_test, np.array([a_test]))[n_test][0]
        pv = np.trapezoid((Tn - Tn_a) * w / (t - a_test), t)
        u_ex = np.pi * chebU_rec(n_test - 1, np.array([a_test]))[n_test - 1][0]
        print(f"  PV-identity test n={n_test}: dense = {pv:.9f}  pi*U_(n-1)(a) = {u_ex:.9f}")

    def dd_solve(load, nup, M=40):
        modes = np.arange(1, M + 1, 2)              # odd Chebyshev modes
        a_col = np.linspace(-0.995, 0.995, 400)
        Umat = chebU_rec(modes.max(), a_col)[modes - 1, :].T
        A = Umat
        rhs = 2 * (1 - nup) * load(a_col * L) * L / G
        coef, *_ = np.linalg.lstsq(A, rhs, rcond=None)
        resid = np.max(np.abs(A @ coef - rhs)) / np.max(np.abs(rhs))
        phiL = np.sum(coef)                         # T_n(1) = 1
        K = phiL * G * np.sqrt(np.pi) / (2 * (1 - nup) * np.sqrt(L))
        # opening du(x) = int_{-L}^{x} b dxi, substitution xi = -L + (x+L) s^2,
        # Gauss-Legendre nodes (no endpoint at the integrable sqrt singularity s=1)
        gls, glw = np.polynomial.legendre.leggauss(2000)
        ss = 0.5 * (gls + 1.0)
        ww = 0.5 * glw
        xq = np.linspace(-L, L, 801)
        du = np.zeros_like(xq)
        for k, x in enumerate(xq):
            xiv = -L + (x + L) * ss * ss
            tv = xiv / L
            Tv = chebT_rec(modes.max(), tv)[modes, :]
            phiv = coef @ Tv
            integrand = phiv * 2 * (x + L) * ss / np.sqrt(
                np.maximum((L - xiv) * (L + xiv), 1e-300))
            du[k] = np.sum(ww * integrand)
        # du(x) = int_{-L}^{x} b dxi = -Delta u(x) with this kernel sign convention
        V = -np.trapezoid(du, xq)
        return K, V, resid

    for tag, load in [("uniform P0=1", lambda x: P0 + 0 * x),
                      ("quadratic 1-(x/L)^2", lambda x: P0 * (1 - (x / L) ** 2))]:
        for nup, lbl in [(nu, "drained nu"), (nu_u, "undrained nu_u")]:
            Kd, Vd, resid = dd_solve(load, nup)
            print(f"DD solve [{tag}, {lbl}]: K = {Kd:.8f}, V = {Vd:.8f}, "
                  f"rel.resid = {resid:.2e}")
    print(f"  expect: uniform load  K = {P0*np.sqrt(np.pi*L):.8f}, "
          f"V(drain) = {np.pi*(1-nu)*P0*L**2/G:.8f}, V(undr) = {np.pi*(1-nu_u)*P0*L**2/G:.8f}")
    print(f"  expect: quadratic load K = {0.5*P0*np.sqrt(np.pi*L):.8f} "
          f"(weight-function integral)")


# ==========================================================================
# CHECK 2: local 1D poroelastic boundary layer (short-time coefficient)
# ==========================================================================
def check2():
    print("\n--- CHECK 2: local 1D poroelastic drainage layer (small-time V growth) ---")
    lhs = 1 - eta * alpha * (1 - 2 * nu) / (1 - nu)
    rhs = (1 - nu_u) / (1 - nu)
    print(f"oedometer identity: 1-eta*alpha(1-2nu)/(1-nu) = {lhs:.12f} ; "
          f"(1-nu_u)/(1-nu) = {rhs:.12f}  -> p_t = c p_yy exactly (diffusivity is c)")
    # FD solve of p_t = c p_yy (c=1, i.e. times in units of 1/c),
    # p(y,0)=2*eta, Dirichlet p(0,t)=0, Neumann dp/dy=0 at y=ymax
    c = 1.0
    ymax = 8.0
    Ny = 2000
    dy = ymax / Ny
    y = np.linspace(0, ymax, Ny + 1)
    p = np.full(Ny + 1, 2 * eta)
    p[0] = 0.0                                        # surface instantly drained
    import scipy.linalg as sla

    def cn_step(p, dt):
        # unknowns q = p[1:] (nodes 1..Ny); node 0 Dirichlet, node Ny mirror (Neumann)
        r = dt * c / dy**2 / 2
        n = Ny
        q = p[1:].copy()
        lower = np.r_[p[0], q[:-1]]
        upper = np.r_[q[1:], q[-2]]                  # mirror at top for last node
        ab = np.zeros((3, n))
        ab[0, 1:] = -r                              # superdiagonal
        ab[1, :] = 1 + 2 * r                         # diagonal
        ab[2, :-1] = -r                              # subdiagonal
        ab[2, -2] = -2 * r                           # last row: mirror entry
        rhs = q + r * (lower - 2 * q + upper)
        qnew = sla.solve_banded((1, 1), ab, rhs)
        out = p.copy()
        out[1:] = qnew
        return out

    targets = list(10.0 ** np.arange(-7, -1.9, 0.5))
    results = []
    t = 0.0
    for tt in targets:
        dt_max = 0.4 * dy * dy / c
        dt = min(dt_max, (tt - t) / 50.0 if tt - t > 0 else dt_max)
        while t < tt * (1 - 1e-14):
            dt = min(dt, tt - t)
            p = cn_step(p, dt)
            t += dt
        integr = np.trapezoid(1.0 - p / (2 * eta), y)
        du_num = (nu_u - nu) / (G * (1 - nu)) * integr
        du_an = 2 * (nu_u - nu) * np.sqrt(c * tt) / (G * np.sqrt(np.pi) * (1 - nu))
        results.append((tt, du_num, du_an))
    print("    t (c*t)    delta u_y(0,t) numeric   analytic 2(nu_u-nu)sqrt(ct)/(G sqrt(pi)(1-nu))  ratio")
    for tt, dn, da in results:
        print(f"  {tt:9.3e}   {dn:16.8e}      {da:16.8e}   {dn/da:8.6f}")
    print("  ratio -> 1 confirms the sqrt(t) surface displacement; then")
    print("  V1(t) = V1_0 + 2L*(2*delta u_y)  =>  A = 8(nu_u-nu)/(sqrt(pi)(1-nu))")
    tt = targets[-1]
    prof = 2 * eta * erf(y / (2 * np.sqrt(c * tt)))
    print(f"  profile check at t={tt:.2e}: max|p_num - 2*eta*erf(y/(2 sqrt(ct)))| = "
          f"{np.max(np.abs(p - prof)):.3e}  (p scale {2*eta:.3f})")


# ==========================================================================
# CHECK 3: loading-2 long-time state
# ==========================================================================
def check3():
    print("\n--- CHECK 3: loading-2 long-time state ---")

    def g_harm(xx, yy):
        z = xx + 1j * yy
        return 1.0 - (z / np.sqrt(z * z - L * L + 0j)).real

    # (a) NEGATIVE check: the bounded FD Laplace box solution (slit p=1, box p=0,
    #     Neumann on y=0 off the slit) vs the decaying harmonic field g(z).
    print(" (a) FD Laplace box vs harmonic candidate g = 1 - Re[z/sqrt(z^2-L^2)]:")
    for X, h in [(6.0, 0.06), (12.0, 0.06)]:
        nx = int(round(2 * X / h)) + 1
        ny = int(round(X / h)) + 1
        x = np.linspace(-X, X, nx)
        y = np.linspace(0, X, ny)
        slit = abs(x) < L - 0.5 * h
        idx = -np.ones((nx, ny), dtype=int)
        u = 0
        for i in range(nx):
            for j in range(ny):
                outer = (i == 0) or (i == nx - 1) or (j == ny - 1)
                if not outer and not (j == 0 and slit[i]):
                    idx[i, j] = u
                    u += 1
        rows, cols, vals, rhs = [], [], [], []
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                if idx[i, j] < 0:
                    continue
                r = len(rhs)
                rows.append(r); cols.append(idx[i, j]); vals.append(-4.0)
                rh = 0.0
                for ii, jj in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                    nb = idx[ii, jj]
                    if nb >= 0:
                        rows.append(r); cols.append(nb); vals.append(1.0)
                    elif jj == 0 and slit[ii]:
                        rh -= 1.0        # slit neighbor below, p=1
                rhs.append(rh)
        for i in range(1, nx - 1):
            if idx[i, 0] >= 0:           # y=0 off the slit: dp/dy = 0
                r = len(rhs)
                rows += [r, r, r]
                cols += [idx[i, 0], idx[i, 1], idx[i, 2]]
                vals += [3.0, -4.0, 1.0]
                rhs.append(0.0)
        A = coo_matrix((vals, (rows, cols)), shape=(len(rhs), u)).tocsr()
        sol = spsolve(A, np.array(rhs))
        P = np.zeros((nx, ny))
        P[idx >= 0] = sol
        P[slit, 0] = 1.0
        for pt in [(0.0, 2.0), (1.5, 0.5)]:
            i = np.argmin(abs(x - pt[0]))
            j = np.argmin(abs(y - pt[1]))
            print(f"   box X={X:4.1f}:  p_FD({pt[0]:.1f},{pt[1]:.1f}) = {P[i, j]: .5f}   "
                  f"g_harmonic = {g_harm(pt[0], pt[1]): .5f}")
        print(f"   box X={X:4.1f}:  min p_FD = {P.min():.4f}  (FD stays in [0,1]; g_harmonic "
              f"is unbounded <0 at the tips)  -> different solutions")
    print("   => no decaying bounded steady state exists; the harmonic candidate is NOT")
    print("      the parabolic limit (also: p_FD at interior points GROWS with box size)")

    # (b) parabolic FD: p_t = lap p with slit held at p=1 -> p -> 1 pointwise
    print(" (b) parabolic FD (slit at p=1, box X=30 held at p=0): fixed points -> 1")
    X, h = 30.0, 0.25
    nx = int(round(2 * X / h)) + 1
    ny = int(round(X / h)) + 1
    x = np.linspace(-X, X, nx)
    y = np.linspace(0, X, ny)
    slit = abs(x) < L - 0.5 * h
    p = np.zeros((nx, ny))
    dt = 0.2 * h * h
    probes = [(0.0, 2.0), (0.0, 5.0), (3.0, 3.0)]
    iprobe = [(np.argmin(abs(x - a)), np.argmin(abs(y - b))) for a, b in probes]
    times = [10.0, 50.0, 200.0]
    out = {T: [] for T in times}
    t = 0.0
    nextT = 0
    nsteps = int(round(times[-1] / dt))
    for step in range(nsteps):
        pn = p.copy()
        pn[1:-1, 1:-1] = p[1:-1, 1:-1] + dt / h**2 * (
            p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4 * p[1:-1, 1:-1])
        # y=0 off the slit: mirror ghost p[i,-1] = p[i,1]
        pn[1:-1, 0] = p[1:-1, 0] + dt / h**2 * (
            p[2:, 0] + p[:-2, 0] + 2 * p[1:-1, 1] - 4 * p[1:-1, 0])
        pn[slit, 0] = 1.0
        p = pn
        t += dt
        while nextT < len(times) and t >= times[nextT] - 1e-9:
            for (ii, jj) in iprobe:
                out[times[nextT]].append(p[ii, jj])
            nextT += 1
    while nextT < len(times):                # flush any remaining probe times
        for (ii, jj) in iprobe:
            out[times[nextT]].append(p[ii, jj])
        nextT += 1
    for T in times:
        vals_ = out[T]
        print(f"   ct = {T:6.1f}:  " + "  ".join(
            f"p({a:.0f},{b:.0f})={v:.4f}" for (a, b), v in zip(probes, vals_)) + "   (-> 1)")
    print("   => p -> 1 pointwise; with sigma -> 0 the final state is uniform swelling")
    print("      u = beta_sw*(x,y): zero relative face displacement -> K2_inf = 0, V2_inf = 0")


# ==========================================================================
def main():
    check1()
    check2()
    check3()

    print("\n" + "=" * 72)
    print("REPORTED VALUES (K in units sqrt(pi L); V in units L^2/G; per unit load)")
    print("=" * 72)
    res = {
        "K1_0": K1_0 / np.sqrt(np.pi * L),
        "K1_inf": K1_inf / np.sqrt(np.pi * L),
        "V1_0": V1_0 * G / L**2,
        "V1_inf": V1_inf * G / L**2,
        "V1_short_coeff": V1_short_coeff,
        "p_faces_0": p_faces_0 / P0,
        "K2_0": K2_0 / np.sqrt(np.pi * L),
        "K2_inf": K2_inf / np.sqrt(np.pi * L),
        "V2_inf": V2_inf * G / L**2,
        "formulas": {
            "K1_0": "K1(0+) = dP*sqrt(pi*L)  (undrained; SIF set by face tractions)",
            "K1_inf": "K1(inf) = dP*sqrt(pi*L)  (drained; same - SIF independent of time)",
            "V1_0": "V1(0+) = pi*(1-nu_u)*dP*L^2/G",
            "V1_inf": "V1(inf) = pi*(1-nu)*dP*L^2/G",
            "V1_short_coeff": "A = 8*(nu_u-nu)/(sqrt(pi)*(1-nu));  V1 ~ V1_0 + A*(L^2/G)*sqrt(c*t/L^2)",
            "p_faces_0": "p_faces(0+) = -eta*(sigma_xx+sigma_yy)|face = 2*(nu_u-nu)/(alpha*(1-2*nu))*dP",
            "K2_0": "0  (nothing has diffused at t=0+; no mechanical load)",
            "K2_inf": "0  (p->1 pointwise, total stress ->0 everywhere)",
            "V2_inf": "V2(inf) = 0  (final state p=1, sigma=0, uniform swelling u=beta_sw*(x,y): "
                      "no relative face displacement; transiently V2<0, closing like ~ -sqrt(ct))",
        },
        "reproduce_command": "python3 solve.py",
    }
    for k in ["K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff", "p_faces_0",
              "K2_0", "K2_inf", "V2_inf"]:
        print(f"  {k:>16s} = {res[k]:.10f}")
    print("  formulas:")
    for k, v in res["formulas"].items():
        print(f"    {k:>16s} : {v}")
    with open("results.json", "w") as f:
        json.dump(res, f, indent=2)
    print("\nwrote results.json")


if __name__ == "__main__":
    main()
