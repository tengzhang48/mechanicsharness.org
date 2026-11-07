"""Grader audit using the copied solver; leaves the original submission intact."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig, solve
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent / "submission"))
from eig_spectral import build_matrices, make_Gfun, mode_value, spectrum
from sim import simulate
from steady import steady_state

root = Path(__file__).resolve().parent.parent
original = json.loads((root / "submission/results.json").read_text())
rerun = json.loads((root / "scratch/submission/results.json").read_text())
print("REPRODUCTION")
for name in ("lambda1", "lambda2", "mode1_shape_half", "K_inf_norm", "relaxation_time_t0"):
    print(name, "original", original[name], "rerun", rerun[name], "difference", rerun[name] - original[name])

print("BASIS CONVERGENCE M=768")
for n in (40, 48, 56, 64, 72):
    lam, vec, _, _, _ = spectrum(n, 768)
    val, c = mode_value(vec, 0, np.array([0.5]))
    print(n, lam[0].real, lam[1].real, val[0] / c[0])

print("QUADRATURE SENSITIVITY N=72")
for m in (384, 1024):
    lam, vec, _, _, _ = spectrum(72, m)
    val, c = mode_value(vec, 0, np.array([0.5]))
    print(m, lam[0].real, lam[1].real, val[0] / c[0])

print("FULL SPECTRUM AND GENERALIZED RESIDUAL N=30 M=384")
a, b, _ = build_matrices(30, 384)
lam, vec = eig(b, a)
order = np.argsort(lam.real)
print("finite", np.isfinite(lam).sum(), "of", len(lam))
print("first eight", lam[order[:8]])
complex_mask = abs(lam.imag) >= 1e-6 * np.maximum(1, abs(lam.real))
print("complex tail", lam[complex_mask])
i = order[0]
res = b @ vec[:, i] - lam[i] * (a @ vec[:, i])
print("absolute residual norm", np.linalg.norm(res))
print("relative residual norm", np.linalg.norm(res) / np.linalg.norm(b @ vec[:, i]))

print("INITIAL-STATE EVOLUTION AND FIRST-MODE AMPLITUDE N=40 M=384")
taus = np.linspace(0, 0.8, 1601)
sigs, z, _ = simulate(40, 384, taus)
mean = np.trapezoid(sigs, z, axis=1)
for tau in (0, 0.02, 0.025, 0.085, 0.1, 0.11, 0.4, 0.5):
    k = np.argmin(abs(taus - tau))
    print("tau", taus[k], "tip", sigs[k, -1], "mean", mean[k],
          "tip*exp(lambda1*tau)", sigs[k, -1] * np.exp(original["lambda1"] * taus[k]))
for ref in (0, 0.025):
    k0 = np.argmin(abs(taus - ref))
    print("reference tau", ref, "mean", mean[k0], "half reached", taus[np.flatnonzero(mean < mean[k0] / 2)[0]],
          "1/e reached", taus[np.flatnonzero(mean < mean[k0] / np.e)[0]])
print("tip < 0.01 reached", taus[np.flatnonzero(sigs[:, -1] < 0.01)[0]])

print("FINITE-TIME TIP INTENSITY FROM MASS CONSERVATION N=40 M=384")
a, b, _ = build_matrices(40, 384)
c0, *_ = np.linalg.lstsq(a, np.ones(a.shape[0]), rcond=None)
generator = -solve(a, b)
w, v = eig(generator)
amplitudes = solve(v, c0)
# Dimensionless wedge density is -integral_0^tau Sigma''' ds.
# Its square-root coefficient is -0.5 * integral G(0) ds;
# K_norm = -pi sqrt(2) times this coefficient.
tip_g = np.array([make_Gfun(j)[0](np.array([0.0]))[0] for j in range(41)])
for tau in (0.02, 0.1, 0.4, 1.0):
    c_int = v @ (np.expm1(w * tau) / w * amplitudes)
    k = np.pi * np.sqrt(2) / 2 * (tip_g @ c_int)
    print("tau", tau, "K_norm", k.real, "imaginary roundoff", k.imag)

print("STEADY CONVERGENCE M=768")
for nc in (32, 48, 64, 96):
    k, *_ = steady_state(nc, 768)
    print(nc, k)

print("POLYNOMIAL EXTRAPOLATION")
data = []
for line in (root / "scratch/submission/eig_poly.stdout.txt").read_text().splitlines()[1:]:
    data.append([float(s) for s in line.split()])
data = np.array(data)
def model(n, limit, coefficient, power):
    return limit + coefficient / n**power
for column in (1, 2):
    fit, _ = curve_fit(model, data[:, 0], data[:, column], p0=(original[f"lambda{column}"], 10, 1), maxfev=10000)
    print("lambda", column, "fit using all printed N", fit)
