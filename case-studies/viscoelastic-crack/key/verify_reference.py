"""Numerical and source checks of the organizer reproduction, not agent grading."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from reference_solver import (
    bulk_dissipation, critical_speed, envelope_kernel_quadrature,
    history_solution, separation, solve_case, strain_convolution, stress,
)


ROOT = Path(__file__).resolve().parent.parent


def main():
    inputs = json.loads((ROOT / "agent/case.json").read_text())
    published = json.loads((ROOT / "private/published_fig9.json").read_text())
    lam, d = inputs["lambda"], inputs["delta_bar"]
    checks = []

    def record(name, error, tolerance, **details):
        checks.append({"name": name, "error": float(error), "tolerance": tolerance,
                       "passed": bool(np.isfinite(error) and error <= tolerance), **details})

    vc = critical_speed(lam, d)
    vc_quad = brentq(
        lambda v: 1 - (1 - lam) * envelope_kernel_quadrature(v) - d,
        0.01, 100.0, xtol=2e-10,
    )
    record("terminal speed series versus direct quadrature", abs(vc - vc_quad), 1e-7,
           series=vc, quadrature=vc_quad,
           scope="Two evaluations of the paper's endpoint equation, not two physical models")
    record("terminal speed versus printed paper", abs(vc - 13.4), 0.05)
    above = solve_case(vc * 1.001, lam, d)
    record("no finite steady root above endpoint", float(above["admissible"]), 0.0)

    for printed in published["cases"]:
        case = solve_case(printed["V"], lam, d)
        v, eps = case["V"], case["epsilon"]
        for key in ("epsilon", "a", "R"):
            record(f"V={v} {key} versus printed Figure 9", abs(case[key]-printed[key]), 0.00051)
        record(f"V={v} cohesive separation equation", abs(separation(eps, v, lam)-d), 1e-9)

        X, Y = np.meshgrid(np.array([-2.0, 0.3, 2.0]), np.pi * np.array([0.1, 0.5, 0.9]), indexing="ij")
        X, Y = X.ravel(), Y.ravel()
        dx = 1e-4
        T = stress(X, Y, eps)
        Tx = (stress(X+dx, Y, eps)-stress(X-dx, Y, eps))/(2*dx)
        Ty = (stress(X, Y+dx, eps)-stress(X, Y-dx, eps))/(2*dx)
        E = strain_convolution(X, Y, case, lam)
        Ex = (strain_convolution(X+dx, Y, case, lam)-strain_convolution(X-dx, Y, case, lam))/(2*dx)
        Ey = (strain_convolution(X, Y+dx, case, lam)-strain_convolution(X, Y-dx, case, lam))/(2*dx)
        record(f"V={v} equilibrium finite difference", np.max(np.abs(Tx.imag+Ty.real)), 2e-6)
        record(f"V={v} strain compatibility finite difference", np.max(np.abs(Ex.real-Ey.imag)), 2e-6)
        record(f"V={v} local SLS constitutive equation", np.max(np.abs(Ex-lam*Tx+(T-E)/v)), 4e-6)

        Xs = np.asarray(inputs["sample_X"])
        Ys = np.pi*np.asarray(inputs["sample_y_over_h"])
        XX, YY = np.meshgrid(Xs, Ys, indexing="ij")
        Ec = strain_convolution(XX, YY, case, lam)
        sol, _, _ = history_solution(Ys, case, lam, rtol=5e-10)
        Eo = lam*stress(XX, YY, eps)+(1-lam)*sol.sol(Xs).T
        record(f"V={v} convolution versus history ODE", np.max(np.abs(Ec-Eo)), 4e-7,
               scope="Different numerical history evaluations sharing the published stress field")

        top = stress(Xs, np.pi, eps)
        top_strain = strain_convolution(Xs, np.pi, case, lam)
        record(f"V={v} fixed-edge tangential stress and strain", max(np.max(np.abs(top.imag)), np.max(np.abs(top_strain.imag))), 1e-10)
        face = stress(np.array([-eps/2, -eps-1]), 1e-9, eps)
        record(f"V={v} cohesive and separated face traction", np.max(np.abs(face.real-np.array([1.0, 0.0]))), 1e-7)
        ahead = stress(32.0, np.pi/2, eps)
        record(f"V={v} relaxed far-ahead stress", abs(ahead-case["a"]), 1e-10)

        energy_runs = [bulk_dissipation(case, lam, d, ny=n) for n in (16, 32, 64)]
        target = case["R"]-1
        error = abs(energy_runs[-1]["value"]-target)/max(1.0, abs(target))
        record(f"V={v} spatial dissipation versus energy release", error, 2e-5,
               target=target, refinements=energy_runs,
               scope="Field integration does not use R-1; both routes share the published model")
        record(f"V={v} energy ordinate refinement", abs(energy_runs[-1]["value"]-energy_runs[-2]["value"])/max(1.0, abs(target)), 2e-5)

    # Domain/tolerance study on the case with the longest wake in the task.
    case = solve_case(10, lam, d)
    ordinary = bulk_dissipation(case, lam, d, ny=64)
    extended = bulk_dissipation(case, lam, d, ny=64, tail=24, rtol=5e-10)
    record("long-wake domain and ODE refinement", abs(extended["value"]-ordinary["value"])/max(1.0, case["R"]-1), 1e-6,
           ordinary=ordinary, extended=extended)

    elastic = [solve_case(v, 1.0, d) for v in inputs["speeds"]]
    record("elastic control R=1 at every speed", max(abs(c["R"]-1) for c in elastic), 1e-10)
    record("elastic control zone independent of speed", max(c["epsilon"] for c in elastic)-min(c["epsilon"] for c in elastic), 1e-10)
    record("elastic control no terminal speed", float(critical_speed(1.0, d) is not None), 0.0)

    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Source reproduction and numerical consistency; no independent reviewer or agent runs",
        "all_passed": all(c["passed"] for c in checks), "checks": checks,
    }
    output = ROOT / "private/verification_results.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"all_passed":result["all_passed"], "checks":len(checks),
                      "failed":[c for c in checks if not c["passed"]]}, indent=2))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
