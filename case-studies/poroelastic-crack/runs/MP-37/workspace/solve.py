"""Reproduce the dimensionless poroelastic crack results and write results.json."""

import argparse
import json
from pathlib import Path

import numpy as np

from factor_check import factor


def calculate():
    nu, nu_u, alpha = 0.2, 0.4, 0.89
    d = nu_u - nu
    a = alpha * (1 - 2 * nu)
    eta = d / (1 - nu_u)
    counts = (65536, 262144, 1048576)
    rows = [factor(eta, count, 0.01) for count in counts]
    # The frequency-zero square-root cusp gives powers of N^(-1/2).
    # With successive fourfold refinement, these weights cancel two terms.
    weights = np.array([1 / 3, -2, 8 / 3])
    n11 = float(weights @ np.array([row["n11"] for row in rows]))
    n21 = float(weights @ np.array([row["n21"] for row in rows]))
    result = {
        "K1_0": n11 / (1 + eta),
        "K1_inf": 1.0,
        "V1_0": float(np.pi * (1 - nu_u)),
        "V1_inf": float(np.pi * (1 - nu)),
        "V1_short_coeff": float(
            4 * d / ((1 + eta) * np.sqrt(np.pi))
            * (4 + np.pi * np.sqrt(2 / eta) * n21)
        ),
        "p_faces_0": 2 * d / a,
        "K2_0": 0.0,
        "K2_inf": -a / (2 * (1 - nu)),
        "V2_inf": float(-np.pi * a / 2),
        "formulas": {
            "normalization": "K_report=K/sqrt(pi*L); V_report=G*V/L^2; unit step loading",
            "definitions": "d=nu_u-nu; a=alpha*(1-2*nu); eta=d/(1-nu_u); n11=N_minus(0;eta)[0,0]; n21=N_minus(0;eta)[1,0]",
            "factor_kernel": "k=abs(xi); m=sqrt(1+xi^2); t=k/(m+k); D=1+eta-2*eta*k*t/m; k_plus=sqrt(-i*xi); k_minus=sqrt(i*xi); m_plus=sqrt(1-i*xi); m_minus=sqrt(1+i*xi); N=[[D,sqrt(2*eta)*t/(m_plus*k_plus)],[sqrt(2*eta)*t/(m_minus*k_minus),1]]; N(0)=diag(1+eta,1)",
            "factor_definition": "N=N_plus*N_minus; N_plus and inverse analytic in Im(xi)>0; N_minus and inverse analytic in Im(xi)<0; N_plus(infinity)=N_minus(infinity)=I; square roots have positive real part",
            "K1_0": "n11/(1+eta) (true-tip r->0 first, then t->0+)",
            "K1_0_bulk_alternative": "1 (t->0+ at fixed distance first, then r->0)",
            "K1_inf": "1",
            "V1_0": "pi*(1-nu_u)",
            "V1_inf": "pi*(1-nu)",
            "V1_short_coeff": "4*d/((1+eta)*sqrt(pi))*(4+pi*sqrt(2/eta)*n21)",
            "p_faces_0": "2*d/a (undrained outer-limit pressure); actual prescribed face trace is 0",
            "K2_0": "0",
            "K2_inf": "-a/(2*(1-nu))",
            "V2_inf": "-pi*a/2",
        },
        "reproduce_command": "OPENBLAS_NUM_THREADS=1 python3 solve.py",
    }
    diagnostics = {
        "factor_grids": rows,
        "extrapolated_n11": n11,
        "extrapolated_n21": n21,
        "extrapolated_norm_error": n11**2 + n21**2 - (1 + eta),
        "linear_extrapolated_K1_0": (2 * rows[-1]["n11"] - rows[-2]["n11"]) / (1 + eta),
        "linear_extrapolated_A": float(
            4 * d / ((1 + eta) * np.sqrt(np.pi))
            * (4 + np.pi * np.sqrt(2 / eta) * (2 * rows[-1]["n21"] - rows[-2]["n21"]))
        ),
        "face_only_A": float(16 * d / ((1 + eta) * np.sqrt(np.pi))),
    }
    return result, diagnostics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", action="store_true")
    args = parser.parse_args()
    results, checks = calculate()
    output = json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n"
    Path("results.json").write_text(output, encoding="ascii")
    print(output, end="")
    if args.diagnostics:
        print(json.dumps(checks, indent=2, sort_keys=True, allow_nan=False))
