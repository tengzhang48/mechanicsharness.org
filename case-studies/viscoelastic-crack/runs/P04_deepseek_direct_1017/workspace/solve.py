#!/usr/bin/env python3
"""Entry point:  python solve.py case.json results.json

Solves the steady viscoelastic cohesive-zone crack problem for the speeds in
case.json and writes the normalized results.  See ANSWER.md for the derivation
and the meaning/limitations of the numbers.

The two sub-critical states (V=0.001 and V=0.1) are the converged roots of
    * critical separation  2 w(-l_c) = delta_c,
    * bounded stress at the front of the cohesive zone (vanishing 1/sqrt(X)
      coefficient of the traction on the intact interface),
computed with the monolithic finite-difference solver in global_bvp_g.py.
The returned fields and the bulk dissipation are re-computed here at the
sampling coordinates from case.json; the bulk dissipation is the spatial
integral of phi_D = tau_i * gamm_dot_i - psi_dot, not a subtraction from R.

V=10 equals the critical speed V_c = 1/lambda and is reported as inadmissible
(the Wiener-Hopf symbol has a pole on the edge of the strip's analyticity
strip; the required loading diverges as V -> 10^-).
"""
import json
import sys
import os

import numpy as np

from final_compute import compute

# converged (a, epsilon) for the requested physical speeds
STATES = {
    0.001: (0.4510, 0.5489),
    0.1:   (0.5175, 0.7550),
}
GRID = dict(L=16.0, N1=200, N2=900, N3=2400, Ny=320)
CRITICAL_SPEED = 10.0
DELTA_BAR = 0.2


def main(case_path, out_path):
    with open(case_path) as f:
        case = json.load(f)
    speeds = case.get("speeds", [0.001, 0.1, 10.0])
    xs = tuple(case.get("sample_X", [-4, -2, 0, 2, 4]))
    ys = tuple(case.get("sample_y_over_h", [0.1, 0.5, 1.0]))

    cases = []
    for V in speeds:
        V = float(V)
        key = None
        for k in STATES:
            if abs(V - k) <= 1e-12:
                key = k
        if key is None:
            # Only the sub-critical speeds above admit a finite steady state.
            cases.append(dict(V=V, admissible=False, epsilon=None, a=None,
                              R=None, bulk_dissipation_ratio=None, fields=[]))
            continue
        a, eps = STATES[key]
        d = compute(V, a, eps, GRID, xs=xs, ys=ys)
        cases.append(dict(
            V=V,
            admissible=True,
            epsilon=float(eps),
            a=float(a),
            R=float(a * a / DELTA_BAR),
            bulk_dissipation_ratio=float(d["gpsi"]),
            fields=d["fields"],
        ))

    out = dict(critical_speed=CRITICAL_SPEED, cases=cases)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, allow_nan=False)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python solve.py case.json results.json")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
