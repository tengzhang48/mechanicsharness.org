#!/usr/bin/env python3
"""
MP-18 -- reproduce command.  Prints lambda1, lambda2, mode1_shape_half,
K_inf_norm, relaxation_time_t0, method.

Model (problem statement, all lengths in h, time in t0 = 4 pi kT h^3/(E* dgb Dgb Omega)):

    sigma_gb(z,t) = sigma0 - (E*/2pi) Int_0^h S(z,zeta) du/dzeta d zeta ,
    S(z,zeta) = 1/(z-zeta) - 1/(z+zeta) - 2 zeta (z-zeta)/(z+zeta)^3   (PV)
    mu = -sigma_gb Omega ,  j = -(dgb Dgb /(kT Omega)) dmu/dz ,
    2 du/dt = -Omega dj/dz .

With s = sigma_gb/sigma0 , x = z/h , tau = t/t0 these reduce to

    d_tau s = (L s)(x) = Int_0^1 S(x,y) s'''(y) dy ,
    s(0,tau) = 0 ,  s'(1,tau) = s''(1,tau) = 0 ,  s(x,0) = 1 .

Method: the PV integral in y is evaluated in the equivalent *potential* form
    Int_0^1 S(x,y) g(y) dy = d_x Int_0^1 K(x,y) g(y) dy ,
    K(x,y) = ln|(x-y)/(x+y)| + 2 x y/(x+y)^2 ,
whose kernel is only *logarithmically* singular; the hat-basis matrix
M_ij = (1/pi) Int_0^1 S(x_i,y) phi_j(y) dy is evaluated with the exact sympy
antiderivative of S (validated to 1e-10 against PV-aware Gauss quadrature).  The
uniform-grid third-derivative D3 carries the boundary conditions; the reduced
operator  A = M D3  is spectrally solved and Richardson-extrapolated in 1/N.
A one-dimensional null space of the discrete D3 (the marginal x-independent
discrete mode) is discarded.
"""
import numpy as np
import stmt_spec
from d3_uniform import build_D3_uniform


def spectrum(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    # physical operator  L = Int_0^1 S(x,y) s'''(y) dy   (no 1/pi: see derivation).
    # stmt_spec.M_stmt carries a factor (1/pi), so multiply it back.
    A = (np.pi * stmt_spec.M_stmt(x, nodes)) @ build_D3_uniform(N)
    ev, V = np.linalg.eig(A)
    order = np.argsort(-ev.real)
    return x, ev[order], V[:, order]


def richardson(vals, Ns):
    Ns = np.array(Ns, float)
    A = np.vstack([np.ones_like(Ns), 1.0 / Ns, 1.0 / Ns ** 2]).T
    return np.linalg.lstsq(A, np.array(vals), rcond=None)[0][0]


def main():
    Ns = [160, 320, 640, 1280]
    lam = {1: [], 2: [], 3: []}
    half = []
    for N in Ns:
        x, ev, V = spectrum(N)
        order = np.argsort(-ev.real)
        nz = [k for k in order if abs(ev[k].real) > 1e-3]     # drop the D3-null mode
        for r in (1, 2, 3):
            lam[r].append(-ev[nz[r - 1]].real)
        v = V[:, nz[0]].real
        v = v / v[-1]                                          # normalise f(1)=1
        half.append(np.interp(0.5, x, v))
        print("N=%4d  lambda1=%.4f  lambda2=%.4f  lambda3=%.4f  f(0.5)=%.5f"
              % (N, lam[1][-1], lam[2][-1], lam[3][-1], half[-1]))

    L1 = richardson(lam[1], Ns)
    L2 = richardson(lam[2], Ns)
    H = richardson(half, Ns)
    print("\nRichardson-extrapolated (1/N, 1/N^2):")
    print("  lambda1 = %.4f" % L1)
    print("  lambda2 = %.4f" % L2)
    print("  f(0.5)  = %.5f   (slowest mode, normalised f(1)=1)" % H)
    print("  [lambda3 at N=%d = %.3f; not extrapolated, N=160 unresolved by 1/N fit]" % (Ns[-1], lam[3][-1]))
    print("  [alternative normalisation  L=(1/pi)Int S s'''  gives lambda1/pi = %.4f]" % (L1 / np.pi))

    # ---- part (b): the t->inf limit is a traction-free crack-like wedge of depth h;
    #      stress intensity factor = edge-crack result (universal, free surface). ----
    K_inf_norm = 1.1215
    # ---- part (c): relaxation time is set by the slowest mode, tau ~ 1/lambda1 in t0 ----
    relax = 1.0 / L1
    print("\nK_inf/(sigma0 sqrt(pi h)) = %.4f  (edge-crack factor of depth h at a free surface)" % K_inf_norm)
    print("relaxation time scale  = 1/lambda1 = %.5f t0" % relax)

    results = {
        "lambda1": round(L1, 4),
        "lambda2": round(L2, 4),
        "mode1_shape_half": round(H, 5),
        "K_inf_norm": round(K_inf_norm, 4),
        "relaxation_time_t0": round(relax, 5),
        "method": ("Galerkin (hat/PV) discretisation of the symmetrized integral operator of the "
                   "statement, d_tau s = Int_0^1 S(x,y) s'''(y) dy, with s(0)=0 and "
                   "s'(1)=s''(1)=0; PV integrals in the potential form d_x Int K s'''; "
                   "spectrum of A=M D3 (uniform FD D3), Richardson-extrapolated in 1/N; "
                   "one discrete D3 null mode discarded."),
        "reproduce_command": "python3 solve_final.py",
    }
    import json
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nwrote results.json")
    return results


if __name__ == "__main__":
    main()