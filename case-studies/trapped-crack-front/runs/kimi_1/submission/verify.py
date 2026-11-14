#!/usr/bin/env python3
"""Independent closed-form checks of results.json."""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

NU = 0.3
R2, R3, Q = 0.6, 0.8, 0.5
BETA = 2.0 * (1.0 - NU) / (2.0 - NU)
D = (1.0 - NU**2) * (1.0 + R2**2) + (1.0 + NU) * R3**2
ETA = 2.0 * R2 * R3 * (BETA * (1.0 + NU) - (1.0 - NU**2)) / D


def closed_form(s):
    """Return g, f, f', |D|f and kappas using real elementary sums."""
    theta = 2.0 * np.pi * np.asarray(s, dtype=float)
    ct, st = np.cos(theta), np.sin(theta)
    den = 1.0 - 2.0 * Q * ct + Q**2
    pc = (Q * ct - Q**2) / den       # sum Q^n cos(n theta)
    ps = Q * st / den                # sum Q^n sin(n theta)
    fc = -0.5 * np.log(den)           # sum Q^n cos(n theta)/n
    fs = np.arctan2(Q * st, 1.0 - Q * ct)  # sum Q^n sin(n theta)/n
    scale = 1.0 / (1.0 + ETA**2)
    f = -scale * (fc - ETA * fs) / np.pi
    fp = 2.0 * scale * (ps + ETA * pc)
    adf = -2.0 * scale * (pc - ETA * ps)
    g = 2.0 * pc
    k1 = -0.5 * adf
    k2 = -0.5 * R2 * adf - R3 * fp
    k3 = -0.5 * R3 * adf + BETA * R2 * fp
    return g, f, fp, adf, (k1, k2, k3)


def main():
    result = json.loads(Path("results.json").read_text())
    b = result["partB"]
    for name in ("f", "kappa1", "kappa2", "kappa3"):
        assert len(b[name]) == 16
    s16 = np.arange(16) / 16.0
    g, f, fp, adf, ks = closed_form(s16)
    for name, values in zip(("f", "kappa1", "kappa2", "kappa3"),
                            (f, *ks)):
        np.testing.assert_allclose(b[name], values, rtol=0.0, atol=2e-12)

    # Verify extrema independently on the symmetric period [-1/2,1/2).
    dense = np.linspace(-0.5, 0.5, 100001)

    def ff(s):
        return closed_form(s)[1]

    def kk3(s):
        return closed_form(s)[4][2]

    for target, reported, objective in (
            ("front_peak_to_peak", b["front_peak_to_peak"], None),
            ("kappa3_max", b["kappa3_max"], kk3)):
        if objective is None:
            values = ff(dense)
            raw = float(np.max(values) - np.min(values))
            np.testing.assert_allclose(reported, raw, rtol=0.0, atol=2e-5)
            imin, imax = np.argmin(values), np.argmax(values)
            h = dense[1] - dense[0]
            rmin = minimize_scalar(ff, bounds=(dense[imin]-h, dense[imin]+h),
                                   method="bounded", options={"xatol": 1e-14})
            rmax = minimize_scalar(lambda s: -ff(s),
                                   bounds=(dense[imax]-h, dense[imax]+h),
                                   method="bounded", options={"xatol": 1e-14})
            refined = -rmax.fun - rmin.fun
            np.testing.assert_allclose(reported, refined, rtol=0.0, atol=2e-12)
        else:
            values = objective(dense)
            imax = int(np.argmax(values))
            h = dense[1] - dense[0]
            opt = minimize_scalar(lambda s: -objective(s),
                                  bounds=(dense[imax]-h, dense[imax]+h),
                                  method="bounded", options={"xatol": 1e-14})
            np.testing.assert_allclose(reported, -opt.fun, rtol=0.0, atol=2e-12)
            np.testing.assert_allclose(b["s_kappa3_max"], opt.x,
                                       rtol=0.0, atol=2e-8)

    # Check algebraic Part A and loading identities independently.
    a = result["partA"]
    np.testing.assert_allclose(a["CG"], -2.0*np.pi, atol=2e-14)
    np.testing.assert_allclose(a["SG"], -2.0*np.pi*ETA, atol=2e-14)
    np.testing.assert_allclose(b["K10_scaled"], np.sqrt((1.0-NU**2)/D),
                               atol=2e-15)
    energy = 2.0*((1.0-NU**2)*(ks[0]+R2*ks[1])
                  +(1.0+NU)*R3*ks[2])/D
    np.testing.assert_allclose(energy, g, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(adf-ETA*fp, -g, rtol=0.0, atol=2e-13)
    period_grid = np.linspace(0.0, 1.0, 100000, endpoint=False)
    assert abs(np.mean(closed_form(period_grid)[1])) < 1e-13
    print("All independent closed-form, equilibrium, energy, JSON-shape,")
    print("loading, and refined-extremum checks passed.")


if __name__ == "__main__":
    main()
