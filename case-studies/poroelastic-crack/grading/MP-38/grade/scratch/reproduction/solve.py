"""Reproducible endpoint formulas and spectral short-time asymptotics."""
import argparse
import json
import math
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy import special, linalg


def spectral(n, cutoff, spacing, qs, nu=0.2, nu_u=0.4):
    k = np.arange(round(cutoff / spacing) + 1) * spacing
    weights = np.ones(len(k))
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    weights *= spacing / (3 * np.pi)
    js = np.empty((2*n, len(k)))
    js[0] = special.j0(k)
    js[1] = special.j1(k)
    low = k < 2*n + 20
    high = ~low
    for order in range(2, 2*n):
        js[order, high] = (2*(order-1)/k[high])*js[order-1, high] - js[order-2, high]
        js[order, low] = special.jv(order, k[low])
    orders = np.arange(1, 2*n, 2)
    phi = np.pi * orders[:, None] * js[1::2] / np.where(k == 0, 1, k)
    phi[:, 0] = 0
    phi[0, 0] = np.pi/2
    psi = np.pi * js[::2]
    signs = (-1.0)**np.arange(n)
    tail_a = np.outer(orders*signs, orders*signs)
    tail_c = np.outer(signs, signs)
    diagonal = np.diag(np.pi*orders/2)
    d = 2*(nu_u-nu)/(1-nu)
    out = []
    for q in qs:
        m = np.sqrt(k*k + q*q)
        aa = -d*k**3/(m*(m+k))
        bb = d*q*q*k/(m*(m+k))
        cc = d*q*q/m
        a = diagonal + (phi*(weights*aa)) @ phi.T
        a -= d*tail_a/(math.sqrt(cutoff**2+q*q)+cutoff)
        b = (phi*(weights*bb)) @ psi.T
        c = (psi*(weights*cc)) @ psi.T
        c += d*q*math.asinh(q/cutoff)*tail_c
        schur = a if d == 0 else a - b @ linalg.solve(c, b.T, assume_a="pos")
        rhs = np.zeros(n)
        rhs[0] = np.pi*(1-nu_u)/2
        u = linalg.solve(schur, rhs, assume_a="pos")
        volume = np.pi*u[0]
        coeff = 2*q*(volume-np.pi*(1-nu_u))/np.sqrt(np.pi)
        tip = np.sum(signs*orders*u)/(1-nu)
        out.append([q, volume, coeff, tip])
    return np.array(out)


def endpoints(short_coeff, initial_tip):
    nu, nu_u, alpha = 0.2, 0.4, 0.89
    eta = alpha*(1-2*nu)/(2*(1-nu))
    return {
        "K1_0": float(initial_tip), "K1_inf": 1.0,
        "V1_0": math.pi*(1-nu_u), "V1_inf": math.pi*(1-nu),
        "V1_short_coeff": float(short_coeff),
        "p_faces_0": 2*(nu_u-nu)/(alpha*(1-2*nu)),
        "K2_0": 0.0, "K2_inf": -eta,
        "V2_inf": -math.pi*alpha*(1-2*nu)/2,
        "formulas": {
            "K1_0": "lim_{q->infinity} sum_{j>=0} (-1)^j*(2*j+1)*u_j(q)/(1-nu); u solves the spectral Schur system specified in ANSWER.md",
            "K1_inf": "1",
            "V1_0": "pi*(1-nu_u)", "V1_inf": "pi*(1-nu)",
            "V1_short_coeff": "(2/sqrt(pi))*lim_{q->infinity} q*(pi*u_0(q)-pi*(1-nu_u)); u solves the spectral Schur system specified in ANSWER.md",
            "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu))",
            "K2_0": "0",
            "K2_inf": "-alpha*(1-2*nu)/(2*(1-nu))",
            "V2_inf": "-pi*alpha*(1-2*nu)/2",
            "spectral_definitions": "q=sqrt(s), m=sqrt(k^2+q^2), d=2*(nu_u-nu)/(1-nu), Phi_j=pi*(2*j+1)*J_(2*j+1)(k)/k, Psi_j=pi*J_(2*j)(k), D_ij=pi*(2*i+1)*delta_ij/2",
            "spectral_A": "A_ij=D_ij-(d/pi)*integral_0^infinity k^3*Phi_i*Phi_j/[m*(m+k)] dk",
            "spectral_B": "B_ij=(d/pi)*integral_0^infinity q^2*k*Phi_i*Psi_j/[m*(m+k)] dk",
            "spectral_C": "C_ij=(d/pi)*integral_0^infinity q^2*Psi_i*Psi_j/m dk",
            "spectral_u": "(A-B*C^(-1)*B^T)*u=(1-nu_u)*(pi/2)*e_0; use infinite-basis limit before q->infinity; at d=0 use A=D without hydraulic block",
            "V1_short": "V1(t)=(L^2/G)*[pi*(1-nu_u)+V1_short_coeff*sqrt(c*t/L^2)+o(sqrt(c*t/L^2))]"
        },
        "reproduce_command": "python3 solve.py",
        "numerical_uncertainty": {"K1_0": 0.001, "V1_short_coeff": 0.01},
        "K1_0_outer_undrained": 1.0
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--explore", action="store_true")
    args = parser.parse_args()
    if args.explore:
        for n in (48, 96, 160):
            values = spectral(n, 40000, 0.25, (100, 200, 400, 800, 1600, 3200))
            print("N", n, "cutoff", 40000, "spacing", 0.25, flush=True)
            print(values, flush=True)
        return
    values = spectral(80, 160000, 0.125, (80, 320, 1280, 5120))
    # Extrapolate the resolved high-frequency response to zero diffusion length.
    x = 1/np.sqrt(values[:, 0])
    coeff = np.polyfit(x, values[:, 2], 3)[-1]
    initial_tip = np.polyfit(x[:3], values[:3, 3], 2)[-1]
    result = endpoints(round(coeff, 2), round(initial_tip, 3))
    with open("results.json", "w", encoding="ascii") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
