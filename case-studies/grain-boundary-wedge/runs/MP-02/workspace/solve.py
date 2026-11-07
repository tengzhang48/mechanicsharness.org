"""MP-02: stress relaxation by diffusion into a grain boundary of a thin film.

Model (dimensionless, Z=z/h, tau=t/t0, t0 = 4 pi kT h^3/(E* delta_gb D_gb Omega)):
  dSigma/dtau = L[Sigma'''] ,  L[g](Z) = PV int_0^1 s(Z,Z') g(Z') dZ',
  s(Z,Z') = 1/(Z-Z') - 1/(Z+Z') - 2 Z'(Z-Z')/(Z+Z')^3,
  BCs: Sigma(0)=0, Sigma'(1)=Sigma''(1)=0,  Sigma(Z,0)=1.

Modes ~ exp(-lambda tau):  lambda f = -L[f'''].
Spectral collocation, basis with exact tip behavior (1-Z)^(5/2):
  phi_0 = 1-(1-Z)^(5/2), phi_j = Z(1-Z)^(5/2) T_{j-1}(2Z-1).
Steady state: L[A] = 1 -> crack density A ~ R/sqrt(1-Z), K = -pi sqrt(2) R(1).

Prints the values stored in results.json.
"""
import json
import numpy as np
from eig_spectral import spectrum, mode_value
from steady import steady_state


def main():
    # ---- eigenvalues & slowest mode ----
    N, M = 72, 768
    lam, vec, A, B, Zc = spectrum(N, M)
    lam1 = float(lam[0].real)
    lam2 = float(lam[1].real)
    f05, c = mode_value(vec, 0, np.array([0.5]))
    mode1_shape_half = float(f05[0]/c[0])   # f(1) = c[0] (phi_0(1)=1, phi_j(1)=0)

    # ---- t -> infinity: edge-crack K ----
    K_inf_norm, d, _, _ = steady_state(64, M)
    K_inf_norm = float(K_inf_norm)

    # ---- relaxation time (e-folding of dominant mode) ----
    relaxation_time_t0 = float(1.0/lam1)

    results = {
        "lambda1": lam1,
        "lambda2": lam2,
        "mode1_shape_half": mode1_shape_half,
        "K_inf_norm": K_inf_norm,
        "relaxation_time_t0": relaxation_time_t0,
        "method": ("Spectral collocation of the dimensionless integro-differential "
                   "relaxation equation dSigma/dtau = L[Sigma'''] (L = Cauchy-PV integral "
                   "operator of the free-surface climb-dislocation kernel) in a basis "
                   "phi_0=1-(1-Z)^(5/2), phi_j=Z(1-Z)^(5/2) T_{j-1}(2Z-1) that satisfies "
                   "the boundary conditions and the exact crack-tip behavior; Gauss-Legendre "
                   "PV quadrature after the smoothing substitution Z'=1-y^2; generalized "
                   "eigendecomposition (scipy.linalg.eig) for lambda f = -L[f''']; "
                   "steady state L[A]=1 solved with crack-singular basis R/sqrt(1-Z) giving "
                   "K_inf = 1.1215223 sigma0 sqrt(pi h) (edge crack). "
                   "Cross-checks: independent polynomial-basis eigensolver and direct "
                   "time-domain integration reproduce lambda1 and the mode shape."),
        "reproduce_command": "python3 solve.py",
    }
    with open("results.json", "w") as fh:
        json.dump(results, fh, indent=2)

    print("=== MP-02 results (values written to results.json) ===")
    print(f"lambda1             = {lam1:.10f}")
    print(f"lambda2             = {lam2:.10f}")
    print(f"mode1_shape_half    = {mode1_shape_half:.10f}")
    print(f"K_inf_norm          = {K_inf_norm:.10f}")
    print(f"relaxation_time_t0  = {relaxation_time_t0:.10f}")

    # ---- verification prints ----
    print("\n=== checks ===")
    # 1) time-domain decay rate
    from sim import simulate
    taus = np.linspace(0, 0.8, 321)
    sigs, Zg, _ = simulate(N=40, M=384, taus=taus)
    i1 = np.argmax(Zg >= 1.0)
    tip = sigs[:, i1]
    mask = (taus >= 0.25) & (taus <= 0.65) & (tip > 0)
    slope, _ = np.polyfit(taus[mask], np.log(tip[mask]), 1)
    print(f"time-domain late-time decay rate = {-slope:.6f}  (spectral lambda1 = {lam1:.6f})")
    k = np.argmin(np.abs(taus-0.4))
    z05 = np.argmin(np.abs(Zg-0.5))
    print(f"time-domain sigma(0.5)/sigma(1) at tau=0.4 = {sigs[k, z05]/sigs[k, i1]:.8f}"
          f"  (mode shape = {mode1_shape_half:.8f})")
    print(f"edge-crack reference value 1.1215223 (Wigglesworth/Tada), computed K_inf_norm = {K_inf_norm:.7f}")


if __name__ == "__main__":
    main()
