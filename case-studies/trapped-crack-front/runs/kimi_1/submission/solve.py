#!/usr/bin/env python3
"""First-order coplanar perturbation of a mixed-mode crack front."""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

NU = 0.3
R2 = 0.6
R3 = 0.8
Q = 0.5
TWO_PI = 2.0 * np.pi

# Rice's coplanar front perturbation, for a Fourier component of the
# advance da(z), is
#   dK1 = -|k| K10 da/2,
#   dK2 = -|k| K20 da/2 - K30 d(da)/dz,
#   dK3 = -|k| K30 da/2 + beta K20 d(da)/dz,
# where beta=2(1-nu)/(2-nu).  Thus the mode-III-to-mode-II term is the
# local rotation, while beta also includes the coplanar elastic coupling.
# Define D and eta so that
#   G0 = K10^2 D/E and dG/G0 = -|D| da + eta d(da)/dz,
# where |D| has symbol |k| and the last term has an epsilon factored out
# when used for Part B.
BETA = 2.0 * (1.0 - NU) / (2.0 - NU)
D = (1.0 - NU**2) * (1.0 + R2**2) + (1.0 + NU) * R3**2
ETA = 2.0 * R2 * R3 * (BETA * (1.0 + NU) - (1.0 - NU**2)) / D


def part_a():
    # da=A cos(kz), with k lambda=2 pi.  All outputs are coefficients of
    # the dimensionless amplitude A/lambda.
    c1 = -np.pi
    s1 = 0.0
    c2 = -np.pi * R2
    s2 = TWO_PI * R3
    c3 = -np.pi * R3
    s3 = -TWO_PI * BETA * R2
    cg = 2.0 * ((1.0 - NU**2) * (c1 + R2 * c2)
                + (1.0 + NU) * R3 * c3) / D
    sg = 2.0 * ((1.0 - NU**2) * (s1 + R2 * s2)
                + (1.0 + NU) * R3 * s3) / D
    return dict(C1=c1, S1=s1, C2=c2, S2=s2,
                C3=c3, S3=s3, CG=cg, SG=sg)


class FrontSeries:
    """Analytic Fourier solution of (|D|-eta d/ds) f=-g."""

    def __init__(self, nmax=96):
        n = np.arange(-nmax, nmax + 1, dtype=int)
        keep = n != 0
        self.n = n[keep]
        an = np.abs(self.n)
        self.gn = Q**an
        self.fn = -self.gn / (TWO_PI * (an - 1j * ETA * self.n))
        # Symbols for derivative and fractional derivative.
        self.dn = TWO_PI * 1j * self.n * self.fn
        self.adn = TWO_PI * an * self.fn

    @staticmethod
    def _sum(n, coeff, s):
        s = np.asarray(s, dtype=float)
        # The terms decay like 2^-|n|; nmax=96 is far beyond roundoff.
        return np.sum(coeff[:, None] * np.exp(1j * TWO_PI * n[:, None]
                                               * s.ravel()[None, :]),
                      axis=0).real.reshape(s.shape)

    def f(self, s):
        return self._sum(self.n, self.fn, s)

    def fp(self, s):
        return self._sum(self.n, self.dn, s)

    def adf(self, s):
        return self._sum(self.n, self.adn, s)

    def g(self, s):
        # Also use the closed form, rather than the Fourier sum, here.
        s = np.asarray(s, dtype=float)
        den = 1.0 - 2.0 * Q * np.cos(TWO_PI * s) + Q**2
        return (1.0 - Q**2) / den - 1.0

    def kappas(self, s):
        adf = self.adf(s)
        fp = self.fp(s)
        return (-0.5 * adf,
                -0.5 * R2 * adf - R3 * fp,
                -0.5 * R3 * adf + BETA * R2 * fp)


def periodic_extrema(func, grid_size=16384):
    """Refine every sampled local extremum on [0,1)."""
    x = np.linspace(0.0, 1.0, grid_size, endpoint=False)
    y = func(x)
    step = 1.0 / grid_size
    # Explicitly refine minima and maxima with the appropriate objective.
    answers = []
    for want_max in (False, True):
        objective = (lambda v: -float(func(v))) if want_max else (lambda v: float(func(v)))
        candidates = []
        for i in range(grid_size):
            yy = y if not want_max else -y
            if yy[i] <= yy[i - 1] and yy[i] <= yy[(i + 1) % grid_size]:
                result = minimize_scalar(objective,
                                         bounds=(x[i] - step, x[i] + step),
                                         method="bounded",
                                         options={"xatol": 1e-14})
                candidates.append((float(result.x % 1.0),
                                   float(func(result.x))))
        answers.append(candidates)
    return answers[0], answers[1]


def main():
    series = FrontSeries()
    s16 = np.arange(16, dtype=float) / 16.0
    k1, k2, k3 = series.kappas(s16)

    minima, maxima = periodic_extrema(series.f)
    fmin = min(minima, key=lambda p: p[1])
    fmax = max(maxima, key=lambda p: p[1])
    k3minima, k3maxima = periodic_extrema(lambda s: series.kappas(s)[2])
    k3max = max(k3maxima, key=lambda p: p[1])
    s_k3max = k3max[0]
    if s_k3max >= 0.5:
        s_k3max -= 1.0

    results = {
        "partA": part_a(),
        "partB": {
            "K10_scaled": float(np.sqrt((1.0 - NU**2) / D)),
            "front_peak_to_peak": float(fmax[1] - fmin[1]),
            "f": series.f(s16).tolist(),
            "kappa1": k1.tolist(),
            "kappa2": k2.tolist(),
            "kappa3": k3.tolist(),
            "kappa3_max": float(k3max[1]),
            "s_kappa3_max": float(s_k3max),
        },
    }
    Path("results.json").write_text(json.dumps(results, indent=2) + "\n")

    # Independent, dense-grid numerical checks.
    s = np.linspace(0.0, 1.0, 20000, endpoint=False)
    f = series.f(s)
    fp = series.fp(s)
    adf = series.adf(s)
    g = series.g(s)
    kap1, kap2, kap3 = series.kappas(s)
    equilibrium_error = np.max(np.abs(adf - ETA * fp + g))
    energy = 2.0 * ((1.0 - NU**2) * (kap1 + R2 * kap2)
                    + (1.0 + NU) * R3 * kap3) / D
    energy_error = np.max(np.abs(energy - g))
    direct_error = np.max(np.abs((-adf + ETA * fp) - g))
    g_fourier = series._sum(series.n, series.gn, s)
    g_error = np.max(np.abs(g_fourier - g))
    print(f"D={D:.15g}, beta={BETA:.15g}, eta={ETA:.15g}")
    print(f"f minimum at s={fmin[0]:.15g}: {fmin[1]:.15g}")
    print(f"f maximum at s={fmax[0]:.15g}: {fmax[1]:.15g}")
    print(f"kappa3 maximum at s={k3max[0]:.15g}: {k3max[1]:.15g}")
    print(f"mean f={np.mean(f):.3e}")
    print(f"max equilibrium residual={equilibrium_error:.3e}")
    print(f"max direct dG/G0 residual={direct_error:.3e}")
    print(f"max reconstructed-energy residual={energy_error:.3e}")
    print(f"max Fourier-g versus closed-form residual={g_error:.3e}")


if __name__ == "__main__":
    main()
