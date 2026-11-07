"""Numerical canonical matrix Wiener-Hopf factorization on a mapped circle."""

import json
import numpy as np


def factor(eta, count, scale=1.0):
    # theta=0 is xi=0; theta=pi represents infinity on the real line.
    theta = 2 * np.pi * np.arange(count) / count
    xi = scale * np.tan(theta / 2)
    k = np.abs(xi)
    m = np.hypot(1, xi)
    t = k / (m + k)
    kernel = np.zeros((count, 2, 2), dtype=complex)
    kernel[:, 0, 0] = 1 + eta - 2 * eta * k * t / m
    kernel[:, 1, 1] = 1
    nonzero = k > 0
    kernel[nonzero, 0, 1] = (
        np.sqrt(2 * eta) * t[nonzero]
        / (np.sqrt(1 - 1j * xi[nonzero]) * np.sqrt(-1j * xi[nonzero]))
    )
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
    n11 = canonical_zero[0, 0].real
    n21 = canonical_zero[1, 0].real
    identity_error = np.max(np.abs(canonical_zero.conj().T @ canonical_zero - kernel[0]))
    residual = np.max(np.abs(project(kernel @ inverse)))
    return {
        "count": count, "scale": scale, "iterations": iteration + 1, "n11": n11, "n21": n21,
        "K1_0": n11 / (1 + eta), "factor_norm_error": float(identity_error),
        "projection_residual": float(residual),
        "imaginary_zero_error": float(np.max(np.abs(canonical_zero.imag))),
    }


if __name__ == "__main__":
    for scale in (0.01, 0.1, 1.0):
        rows = []
        for count in (16384, 65536, 262144):
            row = factor(1 / 3, count, scale)
            rows.append(row)
            print(json.dumps(row, sort_keys=True))
        n11 = 2 * rows[-1]["n11"] - rows[-2]["n11"]
        n21 = 2 * rows[-1]["n21"] - rows[-2]["n21"]
        print(json.dumps({
            "scale": scale, "extrapolated_K1_0": n11 / (4 / 3),
            "extrapolated_A": .6 / np.sqrt(np.pi) * (4 + np.pi * np.sqrt(6) * n21),
            "extrapolated_norm_error": n11**2 + n21**2 - 4 / 3,
        }, sort_keys=True))
