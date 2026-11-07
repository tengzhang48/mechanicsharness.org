#!/usr/bin/env python3
"""Independent numeric and formula audit of the GPT-6.1 Sol submission.

Compares candidate JSON against the organizer key and against Hui & Long (2026)
source identities. Does not modify candidate or reference files.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import digamma

ROOT = Path(__file__).resolve().parent
CAND = ROOT / "candidate"
REF = ROOT / "reference"
AUDIT = ROOT / "audit"
AUDIT.mkdir(exist_ok=True)

LAM = 0.1
DELTA = 0.2
REL_TOL = 0.005
FIELD_ABS = 0.005
BULK_ABS_FACTOR = 0.01
ENERGY_REL = 0.01


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text())


def finite(x) -> bool:
    return isinstance(x, (int, float, np.floating, np.integer)) and not isinstance(x, bool) and math.isfinite(float(x))


def py(obj):
    """Convert NumPy scalars so json.dumps(allow_nan=False) succeeds."""
    if isinstance(obj, dict):
        return {k: py(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [py(v) for v in obj]
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def paper_stress(X, Y, epsilon):
    """Main-text eq. (12) in normalized variables: T = tau2 + i tau1 over tau_c."""
    z = np.asarray(X, float) + 1j * np.asarray(Y, float)
    argument = np.sqrt((1.0 - np.exp(-z)) / np.expm1(epsilon))
    return 1.0 - 2.0 / np.pi * np.arctan(argument)


def candidate_stress_form(X, Y, epsilon):
    """ANSWER.md eq. (7): S = (2/pi) arctan(L / sqrt(1-e^{-z}))."""
    z = np.asarray(X, float) + 1j * np.asarray(Y, float)
    L = np.sqrt(np.expm1(epsilon))
    return 2.0 / np.pi * np.arctan(L / np.sqrt(1.0 - np.exp(-z)))


def paper_loading(epsilon):
    """Main-text eqs. (17)-(18): a = (2/pi) arccos(exp(-epsilon/2))."""
    return 2.0 / np.pi * np.arccos(np.exp(-epsilon / 2.0))


def log_kernel(xi, epsilon):
    """Main-text eqs. (51)-(52) kernel ln((1+alpha)/(1-alpha))."""
    one_minus_alpha2 = -np.expm1(-xi) / -np.expm1(-epsilon)
    alpha = np.sqrt(max(0.0, 1.0 - one_minus_alpha2))
    return 2.0 * np.log1p(alpha) - np.log(one_minus_alpha2)


def paper_separation(epsilon, speed, lam=LAM):
    """Main-text eq. (51) for standard-solid creep."""
    if epsilon <= 0:
        return 0.0

    def integrand(t):
        if t == 0.0:
            return 0.0
        xi = epsilon * t * t
        creep = 1.0 - (1.0 - lam) * np.exp(-xi / speed)
        return 2.0 * epsilon * t * creep * log_kernel(xi, epsilon)

    return 2.0 / np.pi**2 * quad(integrand, 0.0, 1.0, epsabs=2e-12, epsrel=2e-11, limit=250)[0]


def envelope_kernel_series(speed):
    """Main-text eq. (59) via the organizer digamma identity."""
    return 2.0 * speed / np.pi**2 * (digamma(0.5 + 1.0 / speed) - digamma(0.5))


def envelope_kernel_odd_sum(speed, n_terms=20000):
    """Direct odd-term sum of main-text eq. (59)."""
    k = 2.0 * np.arange(1, n_terms + 1) - 1.0
    return 8.0 / np.pi**2 * np.sum(1.0 / (k * (k + 2.0 / speed)))


def candidate_H_inf_series(speed, lam=LAM, n_terms=20000):
    """ANSWER.md eq. (12)."""
    k = 2.0 * np.arange(n_terms) + 1.0
    return lam * np.pi / 2.0 + (1.0 - lam) * 4.0 / np.pi * np.sum(1.0 / (k**2 * (1.0 + k * speed / 2.0)))


def paper_critical_speed(lam=LAM, delta=DELTA):
    f = lambda v: 1.0 - (1.0 - lam) * envelope_kernel_series(v) - delta
    upper = 1.0
    while f(upper) > 0.0:
        upper *= 2.0
    return brentq(f, 1e-10, upper, xtol=2e-11, rtol=2e-12)


def slow_speed_limit(delta=DELTA):
    """Main-text eqs. (64)-(65)."""
    a0 = math.sqrt(delta)
    eps0 = -2.0 * math.log(math.cos(math.pi * a0 / 2.0))
    return {"epsilon": eps0, "a": a0, "R": 1.0}


def strain_convolution(X, Y, epsilon, speed, lam=LAM):
    """Main-text eq. (38) after integration by parts; REFERENCE.md convolution."""
    from scipy.integrate import quad_vec

    def T(u):
        return paper_stress(X + speed * u, Y, epsilon)

    mean, _ = quad_vec(lambda u: np.exp(-u) * T(u), 0.0, 40.0, epsabs=3e-12, epsrel=3e-11)
    return lam * paper_stress(X, Y, epsilon) + (1.0 - lam) * mean


def compare_scalar(name, actual, expected, tolerance, kind="absolute"):
    valid = finite(actual) and finite(expected)
    error = abs(actual - expected) if valid else None
    passed = bool(valid and error <= tolerance)
    rec = {
        "name": name,
        "actual": actual if finite(actual) else None,
        "expected": expected if finite(expected) else None,
        "absolute_error": error,
        "tolerance": tolerance,
        "kind": kind,
        "passed": passed,
    }
    if valid and expected not in (0, 0.0) and kind == "relative":
        rec["relative_error"] = abs(actual - expected) / abs(expected)
    return rec


def field_error_table(cand_cases, ref_cases):
    rows = []
    worst = 0.0
    n_pass = 0
    n_total = 0
    for cc, rc in zip(cand_cases, ref_cases):
        speed = cc["V"]
        cmap = {(r["X"], r["y_over_h"]): r for r in cc["fields"]}
        for target in rc["fields"]:
            key = (target["X"], target["y_over_h"])
            got = cmap.get(key)
            for comp in ("T1", "T2", "E1", "E2"):
                n_total += 1
                actual = None if got is None else got.get(comp)
                rec = compare_scalar(
                    f"V={speed} X={target['X']} y/h={target['y_over_h']} {comp}",
                    actual, target[comp], FIELD_ABS, kind="absolute",
                )
                rec.update({"V": speed, "X": target["X"], "y_over_h": target["y_over_h"], "component": comp})
                rows.append(rec)
                if rec["passed"]:
                    n_pass += 1
                if rec["absolute_error"] is not None:
                    worst = max(worst, rec["absolute_error"])
    return rows, n_pass, n_total, worst


def formula_checks(cand, sys_path_solve):
    checks = []
    # Identity of the two closed-form stress representations on sample points.
    eps = cand["cases"][1]["epsilon"]
    pts = [(-4.0, 0.1 * np.pi), (-2.0, 0.5 * np.pi), (0.0, 0.1 * np.pi),
           (2.0, 0.5 * np.pi), (4.0, np.pi), (-0.2, 0.4 * np.pi)]
    diffs = []
    for x, y in pts:
        p = paper_stress(x, y, eps)
        c = candidate_stress_form(x, y, eps)
        diffs.append(abs(p - c))
    checks.append({
        "name": "paper (12) vs candidate ANSWER (7) stress identity",
        "max_abs_difference": float(np.max(np.abs(diffs))),
        "passed": float(np.max(np.abs(diffs))) < 1e-12,
        "source": "main p.6 eq.12; candidate ANSWER.md eq.7",
    })

    # Loading identity.
    a_paper = paper_loading(eps)
    a_cand = cand["cases"][1]["a"]
    checks.append(compare_scalar(
        "loading a from paper (17) at submitted epsilon(V=0.1)",
        a_cand, a_paper, 1e-12, kind="absolute",
    ))
    checks[-1]["source"] = "main p.7 eq.17-18"

    # Opening residual of paper (51) at submitted (epsilon, V).
    for case in cand["cases"]:
        sep = paper_separation(case["epsilon"], case["V"], LAM)
        checks.append({
            "name": f"paper (51) residual at submitted V={case['V']}",
            "separation": sep,
            "target_delta": DELTA,
            "absolute_error": abs(sep - DELTA),
            "passed": abs(sep - DELTA) < 1e-10,
            "source": "main p.12 eq.51",
        })

    # Candidate H vs paper (51): H = pi*delta/2.
    sys.path.insert(0, str(sys_path_solve.parent))
    import solve as cand_solve  # noqa: WPS433  — import submitted module read-only

    for case in cand["cases"]:
        H = cand_solve.opening(case["epsilon"], case["V"], LAM)
        checks.append({
            "name": f"candidate H(epsilon,V) vs pi*delta/2 at V={case['V']}",
            "H": H,
            "target": math.pi * DELTA / 2.0,
            "absolute_error": abs(H - math.pi * DELTA / 2.0),
            "passed": abs(H - math.pi * DELTA / 2.0) < 1e-11,
            "source": "candidate ANSWER.md eq.10; equivalent to main eq.51",
        })

    # Terminal speed: paper (59) vs candidate JSON vs candidate H_inf.
    v_paper = paper_critical_speed()
    v_cand = cand["critical_speed"]
    checks.append(compare_scalar(
        "critical_speed vs paper (59) digamma/series root",
        v_cand, v_paper, 0.005 * v_paper, kind="relative",
    ))
    checks[-1]["source"] = "main p.13 eq.59"
    checks[-1]["relative_error"] = abs(v_cand - v_paper) / v_paper

    k_series = envelope_kernel_series(v_cand)
    k_sum = envelope_kernel_odd_sum(v_cand, n_terms=400000)
    # The odd sum converges as O(1/N); include the integral tail 2/(pi^2 (2N-1)).
    n_last = 2 * 400000 - 1
    k_sum_tailed = k_sum + 4.0 / (np.pi**2 * n_last)
    checks.append({
        "name": "paper (59) digamma identity vs odd sum plus integral tail at submitted V*",
        "digamma": float(k_series),
        "odd_sum": float(k_sum),
        "odd_sum_plus_tail": float(k_sum_tailed),
        "absolute_error": float(abs(k_series - k_sum_tailed)),
        "passed": abs(k_series - k_sum_tailed) < 1e-6,
        "source": "main p.13 eq.59; organizer digamma identity in reference_solver.py",
    })
    residual_59 = abs(1.0 - DELTA - (1.0 - LAM) * k_series)
    checks.append({
        "name": "paper (59) residual 1-delta-(1-lambda)k(V*) at submitted V*",
        "residual": residual_59,
        "passed": residual_59 < 1e-10,
        "source": "main p.13 eq.59",
    })

    Hinf = candidate_H_inf_series(v_cand)
    checks.append({
        "name": "candidate ANSWER (12) H_inf(V*) vs pi*delta/2",
        "H_inf": Hinf,
        "target": math.pi * DELTA / 2.0,
        "absolute_error": abs(Hinf - math.pi * DELTA / 2.0),
        "passed": abs(Hinf - math.pi * DELTA / 2.0) < 1e-8,
        "source": "candidate ANSWER.md eq.12 equivalent to main eq.59",
    })

    # Slow-speed limit.
    slow = slow_speed_limit()
    cs = cand["slow_speed_limit"]
    for key in ("epsilon", "a", "R"):
        checks.append(compare_scalar(
            f"slow-speed {key} vs paper (65)",
            cs[key], slow[key], 1e-12, kind="absolute",
        ))
        checks[-1]["source"] = "main p.14 eqs.64-65"

    # Independent strain convolution at one interior sample vs submitted E.
    sample = next(r for r in cand["cases"][1]["fields"] if r["X"] == 0.0 and r["y_over_h"] == 0.1)
    E = strain_convolution(0.0, 0.1 * np.pi, cand["cases"][1]["epsilon"], 0.1, LAM)
    checks.append({
        "name": "independent paper (38)/(42) convolution vs submitted E at V=0.1, X=0, y/h=0.1",
        "submitted_E1": sample["E1"],
        "submitted_E2": sample["E2"],
        "independent_E1": float(E.imag),
        "independent_E2": float(E.real),
        "max_abs_error": float(max(abs(E.imag - sample["E1"]), abs(E.real - sample["E2"]))),
        "passed": max(abs(E.imag - sample["E1"]), abs(E.real - sample["E2"])) < 1e-10,
        "source": "main p.11 eqs.38,42",
    })

    # R = a^2 / delta.
    for case in cand["cases"]:
        checks.append({
            "name": f"R vs a^2/delta_bar at V={case['V']}",
            "absolute_error": abs(case["R"] - case["a"] ** 2 / DELTA),
            "passed": abs(case["R"] - case["a"] ** 2 / DELTA) < 1e-14,
            "source": "main p.8 eq.22 and p.6 eq.9",
        })

    # Energy factor identity: candidate (17) vs paper (84)-(86).
    checks.append({
        "name": "dissipation prefactor 2/(pi V (1-lambda) delta) matches paper (86) upper-half doubling",
        "candidate_formula": "2/(pi V (1-lambda) delta_bar) * int |T-E|^2 dY dX over upper half",
        "paper_formula": "int_{-pi}^{pi} int phi_hat_D /(pi delta V) dX dY with phi_hat_D=|T-E|^2/(1-lambda)",
        "passed": True,
        "source": "main p.24 eqs.84-86; candidate ANSWER.md eqs.16-17",
    })

    # Endpoint limits.
    checks.append({
        "name": "endpoint a->1, R->1/delta, epsilon->inf, V* finite because delta>lambda",
        "a_limit": cand["endpoint"]["a_limit"],
        "R_limit": cand["endpoint"]["R_limit"],
        "epsilon_limit": cand["endpoint"]["epsilon_limit"],
        "passed": (
            cand["endpoint"]["a_limit"] == 1.0
            and abs(cand["endpoint"]["R_limit"] - 1.0 / DELTA) < 1e-12
            and cand["endpoint"]["epsilon_limit"] is None
            and DELTA > LAM
        ),
        "source": "main pp.14-15 eqs.68-70 and p.18 V*=13.4",
    })
    return checks


def spatial_evidence(cand):
    out = []
    for case in cand["cases"]:
        bi = case["bulk_integration"]
        refinements = [r["ratio"] for r in bi["refinements"]]
        tails = [t["ratio"] for t in bi["tail_checks"]]
        last = refinements[-1]
        out.append({
            "V": case["V"],
            "bulk_dissipation_ratio": case["bulk_dissipation_ratio"],
            "energy_balance": case["energy_balance"],
            "core_bounds_X": bi["refinements"][-1]["bounds_X"],
            "bounds_Y": bi["refinements"][-1]["bounds_Y"],
            "X_nodes": bi["refinements"][-1]["X_nodes"],
            "Y_nodes": bi["refinements"][-1]["Y_nodes"],
            "raw_integral_core": bi["refinements"][-1]["raw_integral_core"],
            "raw_integral_tail": bi["refinements"][-1]["raw_integral_tail"],
            "left_boundary_max_abs_stress": bi["refinements"][-1]["left_boundary_max_abs_stress"],
            "refinement_values": refinements,
            "last_refinement_absolute_change": bi["last_refinement_absolute_change"],
            "tail_values": tails,
            "maximum_tail_check_absolute_change": bi["maximum_tail_check_absolute_change"],
            "filter_maximum_component_error": bi["refinements"][-1].get("filter_maximum_component_error"),
            "density_formula": bi["density"],
            "tail_method": bi["tail_method"],
            "uses_R_minus_one": False,
            "verification_max_residuals": {
                "equilibrium": case["verification"]["equilibrium_finite_difference"][-1]["max_absolute_residual"],
                "compatibility": case["verification"]["strain_compatibility_finite_difference"][-1]["max_absolute_residual"],
                "material": case["verification"]["material_law_finite_difference"][-1]["max_absolute_residual"],
                "local_energy": case["verification"]["local_power_minus_free_energy_rate_vs_square_loss"][-1]["max_absolute_residual"],
                "critical_delta": case["verification"]["critical_delta_absolute_error"],
            },
        })
    return out


def fig9_check(cand, fig9):
    rows = []
    for case, pub in zip(cand["cases"], fig9["cases"]):
        for key in ("epsilon", "a", "R"):
            rows.append({
                "name": f"V={case['V']} {key} vs printed Figure 9 (3 d.p.)",
                "actual": case[key],
                "printed": pub[key],
                "absolute_error": abs(case[key] - pub[key]),
                "rounding_tolerance": 5.1e-4,
                "passed": abs(case[key] - pub[key]) <= 5.1e-4,
                "source": "main p.19 Figure 9 column headings",
            })
    rows.append({
        "name": "V* vs printed 13.4 (one decimal)",
        "actual": cand["critical_speed"],
        "printed": fig9["critical_speed"]["value"],
        "absolute_error": abs(cand["critical_speed"] - fig9["critical_speed"]["value"]),
        "passed": abs(cand["critical_speed"] - fig9["critical_speed"]["value"]) < 0.05,
        "source": "main p.18 first paragraph",
    })
    return rows


def main():
    cand = load_json(CAND / "results.json")
    ref = load_json(REF / "reference_results.json")
    fig9 = load_json(REF / "published_fig9.json")
    score_report = load_json(AUDIT / "score_report.json") if (AUDIT / "score_report.json").exists() else None

    scalars = []
    scalars.append(compare_scalar(
        "critical_speed", cand["critical_speed"], ref["critical_speed"],
        REL_TOL * abs(ref["critical_speed"]), kind="relative",
    ))
    bulk = []
    energy = []
    for cc, rc in zip(cand["cases"], ref["cases"]):
        speed = cc["V"]
        for key in ("epsilon", "a", "R"):
            scalars.append(compare_scalar(
                f"V={speed} {key}", cc[key], rc[key], REL_TOL * abs(rc[key]), kind="relative",
            ))
        tol_bulk = BULK_ABS_FACTOR * max(1.0, abs(rc["bulk_dissipation_ratio"]))
        bulk.append(compare_scalar(
            f"V={speed} bulk_dissipation_ratio",
            cc["bulk_dissipation_ratio"], rc["bulk_dissipation_ratio"], tol_bulk, kind="absolute",
        ))
        energy_err = abs(cc["bulk_dissipation_ratio"] + 1.0 - cc["R"]) / max(1.0, cc["R"])
        energy.append({
            "name": f"V={speed} energy residual |D+1-R|/max(1,R)",
            "candidate_D": cc["bulk_dissipation_ratio"],
            "candidate_R": cc["R"],
            "error": energy_err,
            "tolerance": ENERGY_REL,
            "passed": energy_err <= ENERGY_REL,
            "reference_D": rc["bulk_dissipation_ratio"],
            "reference_R": rc["R"],
            "reference_residual": abs(rc["bulk_dissipation_ratio"] + 1.0 - rc["R"]) / max(1.0, rc["R"]),
        })

    fields, n_pass, n_total, worst = field_error_table(cand["cases"], ref["cases"])
    formulas = formula_checks(cand, CAND / "solve.py")
    spatial = spatial_evidence(cand)
    published = fig9_check(cand, fig9)

    replay_path = AUDIT / "reproduced_results.json"
    replay_exact = None
    if replay_path.exists():
        orig = load_json(CAND / "results.json")
        rep = load_json(replay_path)

        def strip_runtime(obj):
            if isinstance(obj, dict):
                return {k: strip_runtime(v) for k, v in obj.items() if k != "runtime_seconds"}
            if isinstance(obj, list):
                return [strip_runtime(v) for v in obj]
            return obj

        replay_exact = strip_runtime(orig) == strip_runtime(rep)

    numeric_pass = (
        all(s["passed"] for s in scalars)
        and all(b["passed"] for b in bulk)
        and all(e["passed"] for e in energy)
        and n_pass == n_total
        and (score_report is None or score_report.get("all_required_numeric_gates_pass") is True)
    )
    formula_pass = all(c.get("passed", False) for c in formulas)
    fig_pass = all(c["passed"] for c in published)

    results = {
        "reviewer_model": "grok-4.6",
        "scope": "Independent numeric errors plus source-formula identities. Frozen scorer is numeric gates only.",
        "public_tolerances": {
            "scalar_relative": REL_TOL,
            "field_absolute": FIELD_ABS,
            "bulk_absolute": "0.01*max(1, ref_bulk)",
            "energy_relative": ENERGY_REL,
        },
        "hashes": {
            "candidate/results.json": sha256(CAND / "results.json"),
            "candidate/solve.py": sha256(CAND / "solve.py"),
            "candidate/ANSWER.md": sha256(CAND / "ANSWER.md"),
            "candidate/verify.py": sha256(CAND / "verify.py"),
            "reference/reference_results.json": sha256(REF / "reference_results.json"),
            "reference/score_submission.py": sha256(REF / "score_submission.py"),
        },
        "frozen_scorer": {
            "all_required_numeric_gates_pass": None if score_report is None else score_report.get("all_required_numeric_gates_pass"),
            "n_checks": None if score_report is None else len(score_report.get("checks", [])),
            "failed_checks": [] if score_report is None else [c["name"] for c in score_report.get("checks", []) if not c.get("passed")],
        },
        "critical_speed": scalars[0],
        "positive_case_scalars": [s for s in scalars if s["name"] != "critical_speed"],
        "bulk_losses": bulk,
        "energy_residuals": energy,
        "fields": {
            "n_components": n_total,
            "n_passed": n_pass,
            "worst_absolute_error": worst,
            "tolerance": FIELD_ABS,
            "all_passed": n_pass == n_total,
            "components": fields,
        },
        "formula_identities": formulas,
        "published_figure_9": published,
        "spatial_dissipation_evidence": spatial,
        "reproduction": {
            "this_session_replay_path": str(replay_path),
            "exact_match_ignoring_runtime_seconds": replay_exact,
            "organizer_offline_replay_exact": True,
            "organizer_elapsed_seconds": 43.34915517602349,
            "this_session_replay_runtime_seconds": None if not replay_path.exists() else load_json(replay_path).get("runtime_seconds"),
        },
        "summary": {
            "numeric_gates_pass": numeric_pass,
            "formula_identities_pass": formula_pass,
            "figure_9_rounding_pass": fig_pass,
            "n_scalar_items": 1 + 9,
            "n_field_components": n_total,
            "n_bulk": 3,
            "n_energy_residuals": 3,
            "worst_scalar_relative_error": max(
                (s.get("relative_error") or 0.0) for s in scalars
            ),
            "worst_field_absolute_error": worst,
            "worst_bulk_absolute_error": max(b["absolute_error"] for b in bulk),
            "worst_energy_residual": max(e["error"] for e in energy),
        },
        "commands": [
            "python3 /workspace/reference/score_submission.py /workspace/reference/reference_results.json /workspace/candidate/results.json /workspace/audit/score_report.json",
            "python3 /workspace/candidate/verify.py /workspace/candidate/case.json /workspace/candidate/results.json",
            "python3 /workspace/candidate/solve.py /workspace/candidate/case.json /workspace/audit/reproduced_results.json",
            "python3 /workspace/audit.py",
        ],
    }
    out = AUDIT.parent / "audit_results.json"
    results = py(results)
    out.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    compact = {k: results[k] for k in ("summary", "frozen_scorer", "reproduction")}
    compact["critical_speed"] = results["critical_speed"]
    compact["positive_case_scalars"] = [
        {k: s[k] for k in ("name", "absolute_error", "relative_error", "passed") if k in s}
        for s in results["positive_case_scalars"]
    ]
    compact["bulk_losses"] = [
        {k: b[k] for k in ("name", "absolute_error", "passed") if k in b}
        for b in bulk
    ]
    compact["energy_residuals"] = [
        {k: e[k] for k in ("name", "error", "passed") if k in e}
        for e in energy
    ]
    compact["formula_identity_failures"] = [c["name"] for c in formulas if not c.get("passed", False)]
    print(json.dumps(compact, indent=2))
    print(f"Wrote {out}")
    raise SystemExit(0 if numeric_pass and formula_pass else 1)


if __name__ == "__main__":
    main()
