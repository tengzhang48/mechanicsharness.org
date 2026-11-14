"""
First-order per-mode solver (|k|=1), final formulation (see analysis notes).

Exact WH results (y-sector): b = c- e^{-X} X^{-1/2}, on-crack traction 0, ahead
traction [1/(4(1-nu))] e^{-|x|}|x|^{-3/2}:  A32 = 1/(sqrt(2pi) D1), A12 = -A32,
r_y = -1/2 exactly (Gao & Rice).

Shear sector: b1-column closed forms (pole residue on-crack; erfc/Gamma ahead);
plus-family remainder solved by weighted least squares with the validated Hadamard
columns. Extraction via sqrt(x)-series fits (traction route) and via the sqrt-edge
coefficients (slip route); the identity A12[phi0+] = -A12[b1] (per kernel entry)
makes the two routes equivalent -- they are cross-checked.

b1-column values:
  on-crack shear: e^{-X} * sqrt(pi) cc nu / sqrt(2) * [[-1, i],[i, 1]]
  ahead: T1^{ij}(|x|) = -(sqrt(pi)/pi)[Phat_ij I_P-part + Qhat_ij I_Q]
    Phat_xx = cc nu, Phat_zz = -cc nu, Phat_xz = i cc nu (1+s),
    Qhat_xx = cc, Qhat_zz = cc(1-nu), Qhat_xz = 0,
    I_P = (pi/sqrt2) e^{x} erfc(sqrt(2x)),  I_Q = e^{-x} Gamma(3/2) x^{-3/2},
    xz-part uses  int (1+s) s^-1/2 e^{-(1+s)x}/(s+2) ds = e^{-x} sqrt(pi/x) - I_P.
  => A32^{b1} = (-Q/2, -Q/2, 0) and A12^{b1} = (Q/2, Q/2, -i cc nu) for (xx, zz, xz).
"""
import numpy as np
from scipy.special import erfc, gamma as G
from kernels_x import BasisFun, apply_kernels_matrix, apply_ahead_matrix

NU = 0.3
MU = 1.0
CC = -MU/(2*(1-NU))
D1 = 2*(1-NU)*np.sqrt(2/np.pi)
D3 = 2*np.sqrt(2/np.pi)
S2PI = np.sqrt(2*np.pi)
DVEC = np.array([D1, D3])

def b1_oncrack_shear(X):
    M0 = (np.sqrt(np.pi)*CC*NU/np.sqrt(2.0))*np.array([[-1, 1j], [1j, 1]])
    return np.exp(-X)*M0

def I_P_fun(x):
    return np.exp(-x)*(np.pi/np.sqrt(2.0))*np.exp(2*x)*erfc(np.sqrt(2*x))

def I_Q_fun(x):
    return np.exp(-x)*G(1.5)*x**-1.5

def b1_ahead(entry, x):
    if entry == 'xx':
        return -(np.sqrt(np.pi)/np.pi)*(CC*NU)*I_P_fun(x) - (np.sqrt(np.pi)/np.pi)*CC*I_Q_fun(x)
    if entry == 'zz':
        return -(np.sqrt(np.pi)/np.pi)*(-CC*NU)*I_P_fun(x) - (np.sqrt(np.pi)/np.pi)*CC*(1-NU)*I_Q_fun(x)
    if entry == 'xz':
        IP1 = np.exp(-x)*np.sqrt(np.pi/x) - I_P_fun(x)
        return -(np.sqrt(np.pi)/np.pi)*(1j*CC*NU)*IP1
    if entry == 'yy':
        return -(np.sqrt(np.pi)/np.pi)*CC*I_Q_fun(x)

def b1_A32_A12(entry):
    if entry == 'xz':
        return 0.0, -1j*CC*NU
    Q = {'xx': CC, 'zz': CC*(1-NU), 'yy': CC}[entry]
    return -Q/2.0, +Q/2.0

def make_plus_basis(ML):
    bs = [BasisFun(+0.5, m) for m in range(ML+1)]
    Xs = np.logspace(-6, np.log10(40), 2000)
    for i, b in enumerate(bs):
        bs[i] = BasisFun(b.alpha, b.m, scale=np.max(np.abs(b.eval(Xs))))
    return bs

def ahead_coeffs_reg(basis, coefs, entry, xs=None):
    """A32, A12 of a plus-family combo: T = A32 u^-3 + A12 u^-1 + A0 + ...
    fit T*u and T*u^3 as polynomials in u = sqrt(x)."""
    if xs is None:
        xs = np.array([1e-6, 4e-6, 1.6e-5, 6.4e-5, 2.56e-4, 1.024e-3, 4.1e-3, 1.6e-2])
    Tv = np.array([apply_ahead_matrix(basis, -x, [entry])[entry] @ coefs for x in xs])
    u = np.sqrt(xs)
    M12 = np.stack([np.ones_like(u), u, u**2, u**3, u**4], axis=1)
    s12 = np.linalg.lstsq(M12, Tv*u, rcond=None)[0]
    M32 = np.stack([np.ones_like(u), u**2, u**3, u**4, u**5], axis=1)
    s32 = np.linalg.lstsq(M32, Tv*u**3, rcond=None)[0]
    resid = np.linalg.norm(M12@s12 - Tv*u)/np.linalg.norm(Tv*u)
    return s32[0], s12[0], resid

def solve_shear(ML=6, Ncol=48, Xmin=1e-4, Xmax=25.0, verbose=True):
    basis = make_plus_basis(ML); nb = len(basis)
    Xj = np.logspace(np.log10(Xmin), np.log10(Xmax), Ncol)
    entries = ['xx', 'xz', 'zz']
    Ap = {e: np.zeros((Ncol, nb), complex) for e in entries}
    for j, X in enumerate(Xj):
        d = apply_kernels_matrix(basis, X, entries)
        for e in entries:
            Ap[e][j, :] = d[e]
    B1on = np.stack([b1_oncrack_shear(X) for X in Xj])
    A = np.vstack([np.hstack([Ap['xx'], Ap['xz']]),
                   np.hstack([Ap['xz'], Ap['zz']])])
    rhs = np.zeros((2*Ncol, 2), complex)
    rhs[:Ncol, 0] = -B1on[:, 0, 0]; rhs[Ncol:, 0] = -B1on[:, 1, 0]
    rhs[:Ncol, 1] = -B1on[:, 0, 1]; rhs[Ncol:, 1] = -B1on[:, 1, 1]
    w = np.linalg.norm(np.hstack([A, rhs]), axis=1); w[w == 0] = 1
    sols = []
    for jc in range(2):
        s = np.linalg.lstsq(A/w[:, None], rhs[:, jc]/w, rcond=None)[0]
        res = np.linalg.norm((A@s - rhs[:, jc])/w)/np.linalg.norm(rhs[:, jc]/w)
        sols.append((s, res))
    A12m = np.zeros((2, 2), complex); A32m = np.zeros((2, 2), complex)
    cpm = np.zeros((2, 2), complex)
    for jc in range(2):
        s, res = sols[jc]
        ax, az = s[:nb], s[nb:]
        A32xx_ax, A12xx_ax, _ = ahead_coeffs_reg(basis, ax, 'xx')
        A32xz_ax, A12xz_ax, _ = ahead_coeffs_reg(basis, ax, 'xz')
        A32xz_az, A12xz_az, _ = ahead_coeffs_reg(basis, az, 'xz')
        A32zz_az, A12zz_az, _ = ahead_coeffs_reg(basis, az, 'zz')
        A32xx, A12xx = b1_A32_A12('xx')
        A32zz, A12zz = b1_A32_A12('zz')
        A32xz, A12xz = b1_A32_A12('xz')
        A32m[0, jc] = A32xx_ax + A32xz_az + (A32xx if jc == 0 else 0.0)
        A12m[0, jc] = A12xx_ax + A12xz_az + (A12xx if jc == 0 else A12xz)
        A32m[1, jc] = A32xz_ax + A32zz_az + (A32zz if jc == 1 else 0.0)
        A12m[1, jc] = A12xz_ax + A12zz_az + (A12zz if jc == 1 else A12xz)
        cpm[0, jc] = ax[0]/basis[0].scale - (1.0 if jc == 0 else 0.0)
        cpm[1, jc] = az[0]/basis[0].scale - (1.0 if jc == 1 else 0.0)
        if verbose:
            print(f"  input {jc}: LS res {res:.2e}")
    Dx = 1.0/(S2PI*A32m[0, 0]); Dz = 1.0/(S2PI*A32m[1, 1])
    R = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            R[i, j] = S2PI*(DVEC[j]/2.0)*A12m[i, j]
    if verbose:
        print(f"  A32m real:\n{np.array_str(A32m.real, precision=8)}")
        print(f"  Dx={Dx:.8f} (ex {D1:.8f}), Dz={Dz:.8f} (ex {D3:.8f})")
        print(f"  A12m:\n{np.array_str(A12m, precision=8)}")
        print(f"  R (global):\n{np.array_str(R, precision=8)}")
    return dict(A32m=A32m, A12m=A12m, R=R, cpm=cpm, sols=sols, basis=basis,
                Dx=Dx, Dz=Dz)

def y_sector_exact():
    A32 = 1.0/(S2PI*D1)
    return dict(A32=A32, A12=-A32, r_tract=-0.5, r_slip=-0.5)

def slip_route_R(cpm):
    """R from sqrt-edge coefficients (per unit input channel j: K_j=1, da = 2/D_j,
    a' = i da): c+^x = D2 dK2loc + a' D3 K3 d_jz ; c+^z = D3 dK3loc - a' D2 K2 d_jx
    dK2g = dK2loc + a' K3 d_jz ; dK3g = dK3loc - a' K2 d_jx ; R[i,j] = (D_j/2) dKi^g"""
    R = np.zeros((2, 2), complex)
    for j in range(2):
        Dj = DVEC[j]
        dK2loc = (cpm[0, j] - 2j*(j == 1))/D1
        dK3loc = (cpm[1, j] + 2j*(j == 0))/D3
        dK2g = dK2loc + (2j/D3)*(j == 1)
        dK3g = dK3loc - (2j/D1)*(j == 0)
        R[0, j] = (Dj/2)*dK2g
        R[1, j] = (Dj/2)*dK3g
    return R

if __name__ == "__main__":
    ye = y_sector_exact()
    print("== y-sector (exact WH) ==")
    print(f"  A32 = {ye['A32']:.10f}  r = {ye['r_tract']:.10f}")
    for ML in [5, 6, 8]:
        print(f"== shear sector, ML={ML} ==")
        sh = solve_shear(ML=ML)
        Rsl = slip_route_R(sh['cpm'])
        print("  R slip-route:")
        print(np.array_str(Rsl, precision=8))
        print("  R traction - R slip:")
        print(np.array_str(sh['R'] - Rsl, precision=2))
