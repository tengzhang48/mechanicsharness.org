"""Time-domain simulation: d c/dtau = -A^{-1} B c  (since dSigma/dtau = L[Sigma''']
and B = -L-matrix). Start from uniform Sigma(Z,0)=1 (projected onto the basis).
Tracks decay of the traction -> late-time rate should equal lambda1.
"""
import numpy as np
from scipy.linalg import eig, solve
from eig_spectral import build_matrices, phi_vals


def simulate(N=40, M=384, taus=np.linspace(0, 0.6, 241)):
    A, B, Zc = build_matrices(N, M)
    Minv = -solve(A, B)                    # Mmat = -A^{-1} B
    c0, *_ = np.linalg.lstsq(A, np.ones(len(Zc)), rcond=None)
    w, V = eig(Minv)
    Vi = np.linalg.inv(V)
    amp = Vi @ c0
    # evaluate sigma on a fine physical grid
    Zg = np.linspace(0.0, 1.0, 201)
    Phi = np.zeros((len(Zg), N+1))
    for j in range(N+1):
        Phi[:, j] = phi_vals(j, Zg)
    # sigma(t) = Phi V diag(exp(w t)) amp
    out = []
    for tau in taus:
        c = V @ (np.exp(w*tau)*amp)
        sig = Phi @ c
        out.append(sig.real)
    return np.array(out), Zg, w


if __name__ == "__main__":
    taus = np.linspace(0, 0.8, 321)
    sigs, Zg, w = simulate(N=40, taus=taus)
    # traction at the tip (should be the max at late times)
    i1 = np.argmax(Zg >= 1.0)
    tip = sigs[:, i1]
    # mean traction magnitude
    mean = np.trapezoid(sigs, Zg, axis=1)
    print(f"{'tau':>8} {'sigma(1)':>12} {'mean sigma':>12}")
    for k in range(0, len(taus), 40):
        print(f"{taus[k]:8.3f} {tip[k]:12.6f} {mean[k]:12.6f}")
    # late-time decay rate from tip traction
    mask = (taus >= 0.25) & (taus <= 0.65) & (tip > 0)
    slope, intercept = np.polyfit(taus[mask], np.log(tip[mask]), 1)
    print(f"late-time decay rate from tip traction: {-slope:.6f}  (lambda1 spectral = 11.742668)")
    # normalized profile at late time
    for tau0 in [0.3, 0.5]:
        k = np.argmin(np.abs(taus-tau0))
        prof = sigs[k]/sigs[k, i1]
        z05 = np.argmin(np.abs(Zg-0.5))
        print(f"tau={tau0}: sigma(0.5)/sigma(1) = {prof[z05]:.8f}   (mode1 spectral = 0.80195650)")
    # relaxation time operational: mean falls to 1/e of its early value
    m0 = mean[10]
    k_e = np.argmax(mean < m0/np.e)
    print(f"mean traction reaches 1/e of its tau~0.03 value at tau = {taus[k_e]:.4f}")
    k_h = np.argmax(mean < m0/2.0)
    print(f"mean traction halves at tau = {taus[k_h]:.4f}")
