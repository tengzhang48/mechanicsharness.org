"""Numeric gates only. A pass here is not a reviewed scientific solution."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def score(reference, submission):
    checks = []
    audit = []   # 26 Oct fix (audit of 25 Oct): reported beside the gates announced in CONTRACT.md, not counted in them

    def check(name, passed, **detail):
        checks.append({"name": name, "passed": bool(passed), **detail})

    def audit_check(name, passed, **detail):
        audit.append({"name": name, "passed": bool(passed), **detail})

    def compare(name, actual, expected, tolerance):
        valid = finite_number(actual)
        error = abs(actual-expected) if valid else None
        check(name, valid and error <= tolerance, actual=actual, expected=expected,
              absolute_error=error, tolerance=tolerance)

    if not isinstance(submission, dict):
        check("top-level object", False)
        return {"all_required_numeric_gates_pass": False, "checks": checks,
                "scientific_review": "pending", "study_pass": None}

    vc = reference["critical_speed"]
    if vc is None:
        check("unbounded speed range", "critical_speed" in submission and submission["critical_speed"] is None)
    else:
        compare("critical speed", submission.get("critical_speed"), vc, 0.005*vc)
    actual_cases = submission.get("cases", [])
    if not isinstance(actual_cases, list):
        actual_cases = []
    check("case count", len(actual_cases) == len(reference["cases"]))
    used = set()
    for expected in reference["cases"]:
        speed = expected["V"]
        candidates = [
            (i, c) for i, c in enumerate(actual_cases)
            if isinstance(c, dict) and finite_number(c.get("V"))
            and math.isclose(c["V"], speed, rel_tol=1e-12, abs_tol=1e-14)
        ]
        unique = len(candidates) == 1 and candidates[0][0] not in used
        check(f"V={speed} unique result", unique)
        if not unique:
            continue
        index, actual = candidates[0]
        used.add(index)
        check(f"V={speed} admissibility", type(actual.get("admissible")) is bool
              and actual["admissible"] == expected["admissible"])
        if not expected["admissible"]:
            continue
        for key in ("epsilon", "a", "R"):
            compare(f"V={speed} {key}", actual.get(key), expected[key], 0.005*abs(expected[key]))
        bulk = actual.get("bulk_dissipation_ratio")
        compare(f"V={speed} bulk dissipation", bulk, expected["bulk_dissipation_ratio"],
                0.01*max(1.0, abs(expected["bulk_dissipation_ratio"])))
        R = actual.get("R")
        if finite_number(bulk) and finite_number(R):
            energy_error = abs(bulk+1-R)/max(1.0, R)
            check(f"V={speed} energy balance", energy_error <= 0.01,
                  error=energy_error, tolerance=0.01)
        else:
            check(f"V={speed} energy balance", False, reason="missing finite R or bulk integral")
        # the announced bulk band 0.01 max(1, ref) passes zero at V = 0.001 (ref 0.0073); a relative band does not
        ref = expected["bulk_dissipation_ratio"]
        tol = max(0.05*abs(ref), 2e-4)
        audit_check(f"V={speed} bulk dissipation, 5% band", finite_number(bulk) and abs(bulk-ref) <= tol,
                    actual=bulk, expected=ref, tolerance=tol)
        if finite_number(bulk) and finite_number(R):
            tol = 0.05*abs(R-1)+2e-4
            audit_check(f"V={speed} energy balance, relative", abs(bulk+1-R) <= tol, error=abs(bulk+1-R), tolerance=tol)
        # the submitter's own numbers against Eqs. 18 and 22 of the source: a = (2/pi) acos(exp(-eps/2)), R = a^2/delta_bar
        eps, a = actual.get("epsilon"), actual.get("a")
        # 26 Oct fix: the contract asks for positive epsilon, a, R; test that first, or exp/acos overflow or leave
        # their domain (epsilon = -0.5 or -2000 raised instead of reporting a failure)
        if finite_number(eps) and finite_number(a) and finite_number(R) and eps > 0 and a > 0 and R > 0:
            a_own = 2/math.pi*math.acos(math.exp(-eps/2))
            audit_check(f"V={speed} own a from own epsilon (Eq. 18)", abs(a-a_own) <= 0.005*a, actual=a, implied=a_own)
            R_own = a**2/reference["delta_bar"]
            audit_check(f"V={speed} own R from own a (Eq. 22)", abs(R-R_own) <= 0.005*R, actual=R, implied=R_own)
        else:
            audit_check(f"V={speed} own a, R relations", False, reason="missing, non-finite or non-positive epsilon, a or R")
        rows = actual.get("fields", [])
        if not isinstance(rows, list):
            rows = []
        check(f"V={speed} field row count", len(rows) == len(expected["fields"]))
        for target in expected["fields"]:
            x, y = target["X"], target["y_over_h"]
            matching = [
                r for r in rows if isinstance(r, dict)
                and finite_number(r.get("X")) and finite_number(r.get("y_over_h"))
                and abs(r["X"]-x) < 1e-12 and abs(r["y_over_h"]-y) < 1e-12
            ]
            label = f"V={speed} X={x} y/h={y}"
            check(label+" unique field sample", len(matching) == 1)
            if len(matching) == 1:
                for key in ("T1", "T2", "E1", "E2"):
                    compare(label+" "+key, matching[0].get(key), target[key], 0.005)
    return {
        "all_required_numeric_gates_pass": all(c["passed"] for c in checks),
        "checks": checks, "scientific_review": "pending", "study_pass": None,
        "audit_checks_2026_10_26_pass": all(c["passed"] for c in audit), "audit_checks_2026_10_26": audit,
        "scope": "JSON accuracy only; derivation, integration provenance, reproducibility and access audit require review",
    }


def load_strict(path):
    def parse_float(value):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Nonfinite JSON number: {value}")
        return number

    return json.loads(
        path.read_text(), parse_float=parse_float,
        parse_constant=lambda s: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON number: {s}")),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("submission", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    try:
        result = score(load_strict(args.reference), load_strict(args.submission))
    except (ValueError, OSError) as error:
        result = {"all_required_numeric_gates_pass": False, "format_error": str(error),
                  "scientific_review": "pending", "study_pass": None}
    args.report.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "checks"}, indent=2))
    raise SystemExit(0 if result["all_required_numeric_gates_pass"] else 1)
