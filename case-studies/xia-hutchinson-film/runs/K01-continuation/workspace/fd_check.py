"""Independent numerical check of G_A and G_B by finite differences.

Solves the full 2D membrane-foundation equations on the symmetry half-domain with
traction-free crack faces and symmetry segments, using ghost-node elimination along
the edges (exact tangential identities u1,12 = -u2,11, u2,12 = -nu u1,11 on
traction-free edges; mirror/odd reflections on symmetry edges). G is extracted from
the Eshelby integral  J(c) = int [w - h sigma_{a1} u_{a,1}] dx2 on vertical lines:
path-independent away from the tip; G = (2/h) [J(ahead) - J(behind)].
This route never uses the closed-form tanh results (only the far-field Dirichlet
states, which are the bounded equilibria fixed by the strip ODE - an input to the
model, not the quantity being verified).

Units: s0 = h = E = ell = 1, nu = 0.3, k = h Ebar/ell^2. Analytic targets:
G_A = (1/Ebar) tanh(H/2),  G_B = (1/Ebar) [2 tanh(H/2) - tanh(H)].
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

NU = 0.3
EBAR = 1.0 / (1.0 - NU**2)
S0 = H_TH = 1.0
ELL = 1.0
K = H_TH * EBAR / ELL**2


def wake_u2(x2, S):
    """Relaxed 1D wake on a strip of width S (cracks at both edges), centered form."""
    return -(S0 * ELL / EBAR) * np.sinh((x2 - S / 2.0) / ELL) / np.cosh(S / (2.0 * ELL))


def solve(case, H, n_per=24, Lfac=6.0):
    W = H / 2.0 if case == 'A' else H
    d = W / n_per
    N2 = n_per
    L = Lfac * max(ELL, H)
    n_half = int(round(L / d))
    L = n_half * d
    N1 = 2 * n_half
    i_tip = n_half
    nN = (N1 + 1) * (N2 + 1)

    def pid(i, j, a):
        return 2 * (i * (N2 + 1) + j) + a

    def u_dirichlet(i, j):
        x2 = j * d
        if i == 0:                       # far behind
            S = H if case == 'A' else H  # behind strips have width H in both cases
            return 0.0, wake_u2(x2, S)
        else:                            # far ahead
            if case == 'A':
                return 0.0, 0.0
            return 0.0, wake_u2(x2, 2.0 * H)

    rows, cols, vals = [], [], []
    rhs = []

    def add(r, c, v):
        rows.append(r); cols.append(c); vals.append(v)

    for i in range(N1 + 1):
        for j in range(N2 + 1):
            r1, r2 = pid(i, j, 0), pid(i, j, 1)
            if i == 0 or i == N1:
                v1, v2 = u_dirichlet(i, j)
                add(r1, r1, 1.0); rhs.append(v1)
                add(r2, r2, 1.0); rhs.append(v2)
                continue
            # edge classification
            if j == 0:
                crack = (case == 'B') or (i < i_tip)
            elif j == N2:
                crack = (case == 'B') and (i < i_tip)
            else:
                crack = None
            if crack is None:
                # interior node: full PDE
                # (1) 2 u1,11 + (1-nu) u1,22 + (1+nu) u2,12 - (2/ell^2) u1 = 0
                add(r1, pid(i + 1, j, 0), 2 / d**2); add(r1, pid(i - 1, j, 0), 2 / d**2)
                add(r1, pid(i, j + 1, 0), (1 - NU) / d**2); add(r1, pid(i, j - 1, 0), (1 - NU) / d**2)
                add(r1, r1, -4 / d**2 - 2 * (1 - NU) / d**2 - 2 / ELL**2)
                for (di, dj, sgn) in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
                    add(r1, pid(i + di, j + dj, 1), sgn * (1 + NU) / (4 * d**2))
                # (2) (1-nu) u2,11 + 2 u2,22 + (1+nu) u1,12 - (2/ell^2) u2 = 0
                add(r2, pid(i + 1, j, 1), (1 - NU) / d**2); add(r2, pid(i - 1, j, 1), (1 - NU) / d**2)
                add(r2, pid(i, j + 1, 1), 2 / d**2); add(r2, pid(i, j - 1, 1), 2 / d**2)
                add(r2, r2, -2 * (1 - NU) / d**2 - 4 / d**2 - 2 / ELL**2)
                for (di, dj, sgn) in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
                    add(r2, pid(i + di, j + dj, 0), sgn * (1 + NU) / (4 * d**2))
                rhs.extend((0.0, 0.0))
            elif crack:
                # traction-free edge; jn = +1 if interior is above (j==0), else -1
                jn = 1 if j == 0 else -1
                jb = j + jn
                # u1-eq: (2-nu-nu^2) u1,11 + (1-nu) u1,22[ghost] - (2/ell^2) u1 = 0
                a0 = (2 - NU - NU**2) / d**2
                add(r1, pid(i + 1, j, 0), a0); add(r1, pid(i - 1, j, 0), a0)
                add(r1, r1, -2 * a0 - 2 * (1 - NU) / d**2 - 2 / ELL**2)
                add(r1, pid(i, jb, 0), 2 * (1 - NU) / d**2)
                # ghost term: (1-nu)/d^2 * jn * 2 d u2,1(i,j)
                c = (1 - NU) / d**2 * jn * 2 * d
                add(r1, pid(i + 1, j, 1), c / (2 * d)); add(r1, pid(i - 1, j, 1), -c / (2 * d))
                # u2-eq: -2 nu u2,11 + 2 u2,22[ghost] - (2/ell^2) u2 = 0
                b0 = -2 * NU / d**2
                add(r2, pid(i + 1, j, 1), b0); add(r2, pid(i - 1, j, 1), b0)
                add(r2, r2, -2 * b0 - 4 / d**2 - 2 / ELL**2)
                add(r2, pid(i, jb, 1), 4 / d**2)
                # ghost term: 2/d^2 * jn * 2 d nu u1,1(i,j)
                c2 = 2 / d**2 * jn * 2 * d * NU
                add(r2, pid(i + 1, j, 0), c2 / (2 * d)); add(r2, pid(i - 1, j, 0), -c2 / (2 * d))
                # prestress in sigma22 = s0 + Ebar(u2,2 + nu u1,1) = 0 contributes the
                # constant 2/d^2 * jn * 2 d * s0/Ebar on the LHS -> move to RHS
                rhs.extend((0.0, -jn * 4 * S0 / (EBAR * d)))
            else:
                # symmetry edge: u2 = 0; u1-eq with mirror ghost, u2,12 = 0
                jn = 1 if j == 0 else -1
                jb = j + jn
                a0 = 2 / d**2
                add(r1, pid(i + 1, j, 0), a0); add(r1, pid(i - 1, j, 0), a0)
                add(r1, r1, -2 * a0 - 2 * (1 - NU) / d**2 - 2 / ELL**2)
                add(r1, pid(i, jb, 0), 2 * (1 - NU) / d**2)
                rhs.append(0.0)
                add(r2, r2, 1.0); rhs.append(0.0)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(2 * nN, 2 * nN))
    sol = spla.spsolve(A, np.array(rhs))
    u1 = sol[0::2].reshape(N1 + 1, N2 + 1)
    u2 = sol[1::2].reshape(N1 + 1, N2 + 1)

    # ---- J-integral on vertical lines --------------------------------------
    def d_dx(f, i, j):
        return (f[i + 1, j] - f[i - 1, j]) / (2 * d)

    def d_dy(f, i, j):
        if j == 0:
            return (-3 * f[i, 0] + 4 * f[i, 1] - f[i, 2]) / (2 * d)
        if j == N2:
            return (3 * f[i, N2] - 4 * f[i, N2 - 1] + f[i, N2 - 2]) / (2 * d)
        return (f[i, j + 1] - f[i, j - 1]) / (2 * d)

    Jlines = []
    for i in range(2, N1 - 1):
        f = np.zeros(N2 + 1)
        for j in range(N2 + 1):
            e11 = d_dx(u1, i, j)
            e22 = d_dy(u2, i, j)
            e12 = 0.5 * (d_dy(u1, i, j) + d_dx(u2, i, j))
            ds11 = EBAR * (e11 + NU * e22)
            ds22 = EBAR * (e22 + NU * e11)
            ds12 = EBAR * (1 - NU) * e12
            w = (H_TH * S0 * (e11 + e22)
                 + 0.5 * H_TH * (ds11 * e11 + ds22 * e22 + 2 * ds12 * e12)
                 + 0.5 * K * (u1[i, j]**2 + u2[i, j]**2))
            flux = H_TH * ((S0 + ds11) * e11 + ds12 * d_dx(u2, i, j))
            f[j] = w - flux
        Jlines.append(np.trapezoid(f, dx=d))
    Jlines = np.array(Jlines)
    ilines = np.arange(2, N1 - 1)

    def window(x_lo, x_hi):
        m = (ilines * d - L >= x_lo) & (ilines * d - L <= x_hi)
        return Jlines[m]

    Jb = window(-0.6 * L, -0.25 * L)
    Ja = window(0.25 * L, 0.6 * L)
    G_fd = 2.0 / H_TH * (Ja.mean() - Jb.mean())
    x = H / (2 * ELL)
    G_ex = (np.tanh(x) if case == 'A' else 2 * np.tanh(x) - np.tanh(2 * x)) / EBAR
    return dict(case=case, H=H, d=d, L=L, G_fd=G_fd, G_exact=G_ex,
                rel_err=(G_fd - G_ex) / G_ex,
                J_path_spread_behind=(Jb.max() - Jb.min()) / abs(Jb.mean()),
                J_path_spread_ahead=(Ja.max() - Ja.min()) / max(abs(Ja.mean()), 1e-30),
                J_behind=Jb.mean(), J_ahead=Ja.mean(),
                W_half_behind_exact=(G_ex - (2 * np.tanh(H / ELL) if case == 'B' else 0.0) / EBAR) * H_TH / 2)


def main():
    print('=== FD CHECKS (J-integral route) ===')
    print(f'params: s0=h=E=ell=1, nu={NU}, Ebar={EBAR:.6f}, k={K:.6f}')
    out = []
    for case, H in (('A', 1.0), ('A', 2.0), ('A', 0.5), ('B', 1.0), ('B', 2.0)):
        r = solve(case, H)
        out.append(r)
        print(f"case {case} H={H:3.1f} d={r['d']:.4f}: G_fd={r['G_fd']:.6f} "
              f"G_exact={r['G_exact']:.6f} rel_err={100*r['rel_err']:+.2f}% "
              f"path_spread(b,a)=({r['J_path_spread_behind']:.2e},{r['J_path_spread_ahead']:.2e})")
    # grid convergence: Case A H=1 at d, d/2
    r1 = solve('A', 1.0, n_per=24)
    r2 = solve('A', 1.0, n_per=48)
    rate = (abs(r1['G_fd'] - r1['G_exact']) / max(abs(r2['G_fd'] - r2['G_exact']), 1e-15))
    print(f"convergence A/H=1: err(d)={abs(r1['rel_err']):.3e} err(d/2)={abs(r2['rel_err']):.3e} "
          f"ratio={rate:.2f}")
    # behind-state flux check: J behind ~ W_half(behind) (Case A H=1)
    print(f"J_behind(A,H=1)={r1['J_behind']:.6f} vs -G*h/2={-r1['G_exact']/2:.6f}")
    ok = all(abs(r['rel_err']) < 0.03 for r in out)
    print('FD CHECKS PASS (rel err < 3%)' if ok else 'FD CHECKS: SOME ABOVE 3%')


if __name__ == '__main__':
    main()
