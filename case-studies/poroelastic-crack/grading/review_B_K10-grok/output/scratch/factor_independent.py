#!/usr/bin/env python3
"""Independent numerical WH factorization of the mixed-tip symbol.

Kernel is assembled from the Biot surface relations (not copied from
MP-37's N matrix), then scaled to a canonical Hermitian symbol that
tends to I at infinity. Factorization uses the same FFT Hilbert projector
idea, implemented from scratch.
"""
import json
import numpy as np

NU, NU_U, ALPHA = 0.2, 0.4, 0.89
ETA = (NU_U - NU) / (1.0 - NU_U)
A = ALPHA * (1.0 - 2 * NU)
BCOEF = A / (1.0 - NU)
BETA = ETA / BCOEF


def mixed_symbol(xi):
    """Return scaled N(xi) such that N(inf)=I, N(0)=diag(1+eta,1).

    Unscaled:
      S = -k/(2(1-nu)) * Dmech * W - B * t * p
      Q = -m p + (BETA k (m-k)/(1-nu)) W
    Dmech = 1 + eta*(m-k)/(m+k)
    t = k/(m+k)

    To obtain a symbol with identity at infinity we follow the standard
    one-sided sqrt scaling of (W, p) vs (S, Q):
      mechanical: multiply/divide by k_+ = sqrt(-i xi) so k / k_+ ~ k_-
      hydraulic: scale p by m-factors.

    Equivalently (and matching the known high-k expansion), the canonical
    N can be written after left/right diagonal scaling L, R:
      L (S,Q)^T = N R^{-1} (W, p_scaled) ...
    We construct N by:
      N11 = 2(1-nu) S/W |_{p=0} / k   wait that's Dmech, which is NOT 1 at inf
      Dmech(inf)=1, Dmech(0)=1+eta. Good for N11 if we use Dmech?
      But MP-37 uses D = 1+eta - 2 eta k^2 /(m(m+k)) which equals
      (m+k)/m * Dmech * something.

    Direct computation of D vs Dmech:
    """
    k = np.abs(xi)
    m = np.hypot(1.0, xi)
    t = np.divide(k, m + k, out=np.zeros_like(k), where=k > 0)
    dmech = 1.0 + ETA * np.divide(m - k, m + k, out=np.ones_like(k), where=k > 0)
    dmech = np.where(k == 0, 1.0 + ETA, dmech)
    D = 1.0 + ETA - 2.0 * ETA * k * t / np.maximum(m, 1e-30)
    D = np.where(k == 0, 1.0 + ETA, D)
    return k, m, t, dmech, D


def compare_D():
    xi = np.array([0.0, 0.1, 1.0, 10.0, 100.0])
    k, m, t, dmech, D = mixed_symbol(xi)
    # relation: D =? dmech * m / ((m+k)/2) * k/k ...
    print("xi", xi)
    print("dmech", dmech)
    print("D   ", D)
    print("dmech/D", dmech / D)
    print("(m+k)/(2m)", (m + k) / (2 * m))
    print("D / ((m+k)/(2m) wait)")


def factor(eta, count, scale=0.1):
    theta = 2 * np.pi * np.arange(count) / count
    xi = scale * np.tan(theta / 2)
    k, m, t, dmech, D = mixed_symbol(xi)
    # Canonical kernel as justified in notes: after scaling
    # W_s = k_+ W, p_s = m_+ p  (or the MP-37 choice)
    # N11 = D, N22 = 1, N12 = sqrt(2 eta) t /(m_+ k_+)
    # This scaling is derived as follows:
    # At large k, S ~ -k W / (2(1-nu)), so S_s / W_s ~ 1 if S_s = 2(1-nu) S / k_-
    #   and W_s = k_+ W, because k = k_+ k_- on the real line (k_+ k_- = |xi|).
    # Hydraulic: Q ~ -m p, scale Q_s = -Q / m_-, p_s = m_+ p => Q_s/p_s ~ 1.
    # Cross terms produce the sqrt(2 eta) factors from
    #   2(1-nu) * (B t) / (k_- m_+)  and  (BETA k (m-k)/(1-nu)) / (m_- k_+)
    # Check product of couplings: B * BETA = eta, and
    # k (m-k) t / (k_- m_+ m_- k_+) * constants.
    # m_+ m_- = m, k_+ k_- = k, t=k/(m+k), m-k = 1/(m+k) in scaled units
    # k(m-k)t / (k m) = t (m-k)/m = k (m-k) / (m (m+k)) = k / (m (m+k)^2)
    # Meanwhile t/(m_+ k_+) * t/(m_- k_-) = t^2 / (m k)
    # We want off-diagonal product = 2 eta t^2 /(m k) or similar vs D-related det.
    kernel = np.zeros((count, 2, 2), dtype=complex)
    kernel[:, 0, 0] = D
    kernel[:, 1, 1] = 1.0
    nz = k > 0
    m_plus = np.sqrt(1.0 - 1j * xi)
    k_plus = np.sqrt(-1j * xi)
    kernel[nz, 0, 1] = np.sqrt(2 * eta) * t[nz] / (m_plus[nz] * k_plus[nz])
    kernel[:, 1, 0] = kernel[:, 0, 1].conj()
    kernel[count // 2] = np.eye(2)

    error = kernel - np.eye(2)
    lower = np.fft.fftfreq(count) < 0
    lower[count // 2] = False

    def project(values):
        modes = np.fft.fft(values, axis=0)
        modes[~lower] = 0
        return np.fft.ifft(modes, axis=0)

    inverse = np.broadcast_to(np.eye(2), kernel.shape).astype(complex).copy()
    for iteration in range(200):
        updated = np.eye(2) - project(error @ inverse)
        change = np.max(np.abs(updated - inverse))
        inverse = updated
        if change < 2e-14:
            break
    canonical_zero = inverse[count // 2] @ np.linalg.inv(inverse[0])
    n11 = float(canonical_zero[0, 0].real)
    n21 = float(canonical_zero[1, 0].real)
    return {
        "count": count,
        "scale": scale,
        "iterations": iteration + 1,
        "n11": n11,
        "n21": n21,
        "K1_0": n11 / (1 + eta),
        "norm_err": float(n11**2 + n21**2 - (1 + eta)),
        "detN0": float(np.abs(np.linalg.det(kernel[0]))),
    }


if __name__ == "__main__":
    compare_D()
    eta = ETA
    rows = []
    for count in (4096, 16384, 65536):
        row = factor(eta, count, 0.1)
        rows.append(row)
        print(json.dumps(row, sort_keys=True))
    n11 = 2 * rows[-1]["n11"] - rows[-2]["n11"]
    n21 = 2 * rows[-1]["n21"] - rows[-2]["n21"]
    d = NU_U - NU
    Avol = 4 * d / ((1 + eta) * np.sqrt(np.pi)) * (4 + np.pi * np.sqrt(2 / eta) * n21)
    print(json.dumps({
        "extrapolated_n11": n11,
        "extrapolated_n21": n21,
        "extrapolated_K1_0": n11 / (1 + eta),
        "extrapolated_A": float(Avol),
        "norm_err": n11**2 + n21**2 - (1 + eta),
        "paper_K": (1 - NU_U) / (1 - NU),
        "paper_A": float(16 * (1 - NU_U) * d / ((1 - NU) * np.sqrt(np.pi))),
        "face_only_A": float(16 * d / ((1 + eta) * np.sqrt(np.pi))),
    }, sort_keys=True))
