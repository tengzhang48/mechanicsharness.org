"""Reproduce all reported numbers for the GB-diffusion wedge eigenproblem (MP-47).

Eigenproblem (dimensionless, z = h x, xi = h xi):
    lambda * phi(x) = - int_0^1 S(x,xi) phi'''(xi) d xi
    S(x,xi) = 1/(x-xi) - 1/(x+xi) - 2 xi (x-xi)/(x+xi)^3   (PV at xi=x)
    BCs: phi(0)=0, phi'(1)=0, phi''(1)=0.

Method: Chebyshev-Lobatto collocation for the *values* of phi; the singular
integral operator is applied through a high-accuracy Gauss-Legendre quadrature
plus synthetic Chebyshev interpolation (op2.build_Op2, verified to ~1e-12),
decoupling quadrature accuracy from collocation.  The three BCs are imposed
through the null space of C=[e_N ; D1[0] ; D2[0]] (node 0 is x=1); one
redundant interior equation (at the node nearest x=1) is dropped to make the
pencil square.  An independent overdetermined SVD scan (scan3.py) locates the
same ladder.

Run:  python3 solve.py
"""
import json
import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD
from op2 import build_Op2

K_INF_NORM = 1.1215          # edge crack in a half-plane under remote tension


def solve_drop(N, ngl=700, drop='start'):
    w0, xc = clencurt(N); x = (xc + 1) / 2.0
    D0, _ = chebD(N); Dz = 2 * D0; D1 = Dz; D2 = Dz @ Dz; D3 = Dz @ Dz @ Dz
    Op = build_Op2(x, N, ngl)          # (N-1) x (N+1), rows = interior nodes 1..N-1
    Aint = -Op @ D3
    eN = np.zeros(N + 1); eN[N] = 1.0  # phi(0)=0 at node N (x=0)
    C = np.vstack([eN, D1[0, :], D2[0, :]])   # node 0 is x=1
    B = null_space(C)
    A = Aint @ B
    Mall = B[1:N, :]
    if drop == 'start':                # drop equation at node nearest x=1
        keep = np.arange(1, N - 1)
    elif drop == 'mid':
        keep = np.concatenate([np.arange(0, (N - 2) // 2),
                               np.arange((N - 2) // 2 + 1, N - 1)])
    lam, vec = eig(A[keep, :], Mall[keep, :])
    return x, B, lam, vec


def first_modes(N, ngl=700, drop='start', nmax=2):
    x, B, lam, vec = solve_drop(N, ngl, drop)
    order = np.argsort(lam.real)
    pos = [k for k in order if lam[k].real > 1e-6][:nmax]
    out = []
    for k in pos:
        v = (B @ vec[:, k]).real
        phi = v / v[0]                       # f(1) = 1
        out.append((lam[k].real, x, phi))
    return out


if __name__ == "__main__":
    # --- convergence of the two slowest rates (drop=start) -------------
    print("convergence check (drop='start'):")
    for N in [200, 250, 300, 350, 400]:
        m = first_modes(N, nmax=2)
        print("  N=%3d  lambda1=%.5f  lambda2=%.4f" % (N, m[0][0], m[1][0]))

    # --- adopted values at N=300 --------------------------------------
    N = 300
    (lam1, x, phi1), (lam2, _, _) = first_modes(N, nmax=2)
    fhalf = phi1[int(np.argmin(np.abs(x - 0.5)))]
    relax = 1.0 / lam1

    print("\nN=300 mode-1 shape, f(1)=1:")
    for xv in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
        i = int(np.argmin(np.abs(x - xv)))
        print("   z/h=%.2f  f=%.6f   (lambda1 contribution)" % (x[i], phi1[i]))

    results = {
        "lambda1": round(float(lam1), 4),
        "lambda2": round(float(lam2), 4),
        "mode1_shape_half": round(float(fhalf), 5),
        "K_inf_norm": K_INF_NORM,
        "relaxation_time_t0": round(float(relax), 4),
        "method": ("Chebyshev-Lobatto collocation of the singular integral "
                   "operator lambda*phi=-int S phi'''; Gauss-Legendre quadrature "
                   "with synthetic Chebyshev interpolation (~1e-12); BCs imposed "
                   "via null space of [phi(0),phi'(1),phi''(1)]; confirmed by an "
                   "overdetermined SVD (sigma_min) scan and by N-convergence."),
        "reproduce_command": "python3 solve.py",
    }
    print("\n" + json.dumps(results, indent=2))
    with open("results.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print("\nwrote results.json")