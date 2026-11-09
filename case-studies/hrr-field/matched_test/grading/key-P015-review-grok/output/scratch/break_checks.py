"""Attempts to break crack_tip_checks: false negatives, false positives, grid/tolerance.
Work only on a copy of the checker. Run: python3 scratch/break_checks.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crack_tip_checks import (
    check_field,
    williams_exponent_symmetric,
    hrr_crack_exponent,
    format_report,
)


def failed(report):
    return {r["check"] for r in report if r["passed"] is False}


def passed(report):
    return {r["check"] for r in report if r["passed"] is True}


def by_name(report):
    return {r["check"]: r for r in report}


def dump(title, report, extra=""):
    fails = sorted(failed(report))
    print("\n===", title, "===")
    if extra:
        print(extra)
    print("FAIL:", fails or "(none)")
    for r in report:
        flag = {True: "PASS", False: "FAIL", None: "n/a"}[r["passed"]]
        print(f"  {flag:4} {r['check']}: value={r['value']} tol={r['tol']}")


def williams_crack(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    k = 1 / np.sqrt(2 * np.pi)
    return (
        th,
        k * (1.25 * np.cos(th / 2) - 0.25 * np.cos(1.5 * th)),
        k * (0.75 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th)),
        k * (0.25 * np.sin(th / 2) + 0.25 * np.sin(1.5 * th)),
    )


def williams_notch(alpha, n_theta=4001, mode="I"):
    if mode == "I":
        p = williams_exponent_symmetric(alpha)
        lam = 1 - p
        s = lam + 1
        th = np.linspace(0, alpha, n_theta)
        B = -np.cos((lam + 1) * alpha) / np.cos((lam - 1) * alpha)
        F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
        Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
        Fpp = -(lam + 1) ** 2 * np.cos((lam + 1) * th) - B * (lam - 1) ** 2 * np.cos(
            (lam - 1) * th
        )
        return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p
    # mode II: F = A sin((lam+1)th) + B sin((lam-1)th)
    # characteristic sin(2 lam a) - lam sin(2a) = 0
    from scipy.optimize import brentq

    def g(lam):
        return np.sin(2 * lam * alpha) - lam * np.sin(2 * alpha)

    grid = np.linspace(1e-6, 1.0, 20001)
    v = g(grid)
    lam = None
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            lam = brentq(g, grid[i], grid[i + 1])
            break
    if lam is None:
        raise ValueError(f"no mode II root for alpha={alpha}")
    p = 1 - lam
    s = lam + 1
    th = np.linspace(0, alpha, n_theta)
    # F(alpha)=0: A sin((lam+1)a)+B sin((lam-1)a)=0
    A = 1.0
    den = np.sin((lam - 1) * alpha)
    B = -A * np.sin((lam + 1) * alpha) / den
    F = A * np.sin((lam + 1) * th) + B * np.sin((lam - 1) * th)
    Fp = A * (lam + 1) * np.cos((lam + 1) * th) + B * (lam - 1) * np.cos((lam - 1) * th)
    Fpp = -A * (lam + 1) ** 2 * np.sin((lam + 1) * th) - B * (lam - 1) ** 2 * np.sin(
        (lam - 1) * th
    )
    return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p


def williams_modeII_crack(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    k = 1 / np.sqrt(2 * np.pi)
    # standard mode II Williams
    srr = k * (-0.25 * np.sin(th / 2) + 0.75 * np.sin(1.5 * th))
    stt = k * (-0.75 * np.sin(th / 2) - 0.25 * np.sin(1.5 * th))
    srt = k * (0.25 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th))
    return th, srr, stt, srt


def main():
    results = []

    # --- A. correct fields that might fail (false positives) ---
    th, a, b, c = williams_crack()
    rep = check_field(th, a, b, c, p=0.5, law="linear", nu=0.3, plane="strain")
    dump("Williams crack plane strain (expect pass)", rep)
    results.append(("FP?", "williams crack PE", failed(rep)))

    rep = check_field(th, a, b, c, p=0.5, law="linear", nu=0.3, plane="stress")
    dump("Williams crack plane stress (expect pass)", rep)
    results.append(("FP?", "williams crack PS", failed(rep)))

    # n=1 power law (incompressible linear) on Williams field — should FAIL
    # because Williams is compressible unless nu=0.5 and the angular functions
    # of n=1 HRR are NOT the same as nu=0.3 Williams.
    rep = check_field(th, a, b, c, p=0.5, n=1, law="power", plane="strain")
    dump("Williams crack declared as n=1 power (expect compatibility fail)", rep)
    results.append(("handled?", "williams as n=1 power", failed(rep)))

    # nu=0.5 linear on Williams — plane strain incompressible linear should match
    # power n=1? Linear nu=0.5 plane strain: eps_zz=0, volumetric?
    # Hooke: szz=nu(srr+stt)=0.5(srr+stt), same as incompressible.
    # But strain: e=(1+nu)s - nu tr with E=1. Power n=1: eps=1.5 s_ij (deviatoric)
    # 1.5 * s_ij vs (1+0.5)s - 0.5 tr = 1.5 s - 0.5 tr
    # For incompressible, tr_stress = 3m, s_ij = sigma - m, 1.5 s vs 1.5(sigma-m)=1.5 sigma - 0.5 tr
    # Same IF hydrostatic part handled. Linear e(s)=(1.5)s - 0.5 tr. Power: 1.5*(s-m)=1.5s-0.5 tr.
    # So strains MATCH for nu=0.5. Williams field is valid for any nu (eigenvalue independent of nu).
    # So Williams + law=linear nu=0.5 AND Williams + law=power n=1 should both pass.
    # Wait: power n=1 uses se**(n-1)=1, eps=1.5 s_ij. Linear nu=0.5 as above matches.
    # Williams field satisfies biharmonic, hence compatibility for linear. So n=1 power on Williams
    # SHOULD pass if plane strain (incompressible).
    # Earlier dump will tell us.

    for npts in [51, 101, 401, 1001]:
        thc, ac, bc, cc = williams_crack(n_theta=npts)
        rep = check_field(thc, ac, bc, cc, p=0.5, law="linear", nu=0.3, plane="strain")
        dump(f"Williams crack coarse grid n_theta={npts}", rep)
        results.append((f"grid {npts}", failed(rep)))

    # other notch angles
    for alpha in [0.6 * np.pi, 0.75 * np.pi, 0.9 * np.pi, np.pi * 0.99]:
        try:
            thn, an, bn, cn, p = williams_notch(alpha)
            rep = check_field(thn, an, bn, cn, p=p, law="linear", alpha=alpha)
            dump(f"Williams notch alpha={alpha:.4f} p={p:.6f}", rep)
            results.append((f"notch {alpha:.3f}", failed(rep)))
        except Exception as e:
            print("EXC notch", alpha, e)

    # mode II crack
    th2, a2, b2, c2 = williams_modeII_crack()
    rep = check_field(th2, a2, b2, c2, p=0.5, law="linear", mode="II", plane="strain")
    dump("Williams mode II crack (displacement/J skipped)", rep)
    results.append(("mode II crack", failed(rep)))

    # mode II notch (150 deg, has a root)
    alpha = np.deg2rad(150)
    try:
        thn, an, bn, cn, p = williams_notch(alpha, mode="II")
        rep = check_field(
            thn, an, bn, cn, p=p, law="linear", alpha=alpha, mode="II"
        )
        dump(f"Williams mode II notch 150deg p={p:.6f}", rep)
        results.append(("mode II notch", failed(rep)))
    except Exception as e:
        print("EXC mode II notch", type(e), e)

    # Williams field with wrong nu still passes equilibrium (nu not in eq)
    rep = check_field(*williams_crack(), p=0.5, law="linear", nu=0.49, plane="strain")
    dump("Williams crack nu=0.49 (should still pass; linear eigen independent of nu)", rep)
    results.append(("nu 0.49", failed(rep)))

    # --- B. wrong fields that might pass (false negatives) ---

    # 1. Higher-order Williams term (r^0 T-stress angular? or lambda=1.5, p=-0.5)
    # Mode I next eigenvalue for a crack: lam=1.5, p=-0.5 (bounded)
    th = np.linspace(-np.pi, np.pi, 4001)
    lam = 1.5
    # F = cos((lam+1)th) + B cos((lam-1)th); for crack alpha=pi,
    # B = -cos((lam+1)pi)/cos((lam-1)pi)
    B = -np.cos((lam + 1) * np.pi) / np.cos((lam - 1) * np.pi)
    F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
    Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
    Fpp = -(lam + 1) ** 2 * np.cos((lam + 1) * th) - B * (lam - 1) ** 2 * np.cos(
        (lam - 1) * th
    )
    s = lam + 1
    p = 1 - lam  # -0.5
    srr, stt, srt = Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp
    rep = check_field(th, srr, stt, srt, p=p, law="linear")
    dump(f"Higher-order Williams lam=1.5 p={p} (admissible but not leading)", rep)
    results.append(("higher-order Williams", failed(rep)))

    # 2. T-stress (uniform sigma_xx): p=0, Srr=cos^2, Stt=sin^2, Srt=-sin cos
    srr = np.cos(th) ** 2
    stt = np.sin(th) ** 2
    srt = -np.sin(th) * np.cos(th)
    rep = check_field(th, srr, stt, srt, p=0.0, law="linear")
    dump("T-stress p=0 (not leading, faces: sigma_theta=sin^2 not 0 at +-pi? at pi sin=0)", rep)
    results.append(("T-stress", failed(rep)))

    # 3. Hydrostatic angular field that satisfies faces? hard.
    # Uniform hydrostatic: Srr=Stt=1, Srt=0, p=0. Faces Stt=1 !=0 FAIL faces.

    # 4. Scale-invariant wrong constitutive: take Williams, declare power n=3 p=0.5
    # equilibrium OK (same p), faces OK, symmetry OK, energy: q=-1.5, q-p+2=0.0 NOT >0?  -1.5-0.5+2=0
    # finite energy: q-p+2 = -np -p +2 = -p(n+1)+2 = -0.5*4+2=0, NOT >0 so FAILS energy.
    rep = check_field(*williams_crack(), p=0.5, n=3, law="power", plane="strain")
    dump("Williams declared as n=3 power p=0.5 (wrong law)", rep)
    results.append(("williams as n=3", failed(rep)))

    # 5. Wrong field but p and n related so energy and J-exponent pass:
    # take Williams (p=0.5) declare n=1 power. Energy: p(n+1)=1<2. J exponent target 1/2.
    # Compatibility of power n=1 vs Williams: as argued, should pass in plane strain.
    # That's actually a CORRECT field for incompressible linear, not a false negative.

    # 6. Linear field with p slightly off but within 1e-3 eq tol while exponent check is 1e-6
    # Angular eq residual for wrong p: the field satisfies (1-p0)Srr - Stt + Srt'=0
    # residual at p: (p0-p)Srr, scale ~ |p-p0|.
    # So eq residual /scale ~ |dp|. tol=1e-3 so dp=5e-4 might pass eq but fail exponent (1e-6).
    th, a, b, c = williams_crack()
    rep = check_field(th, a, b, c, p=0.5 + 5e-4, law="linear")
    dump("Williams with p=0.5005 (eq maybe pass, exponent should fail)", rep)
    results.append(("p off 5e-4", failed(rep)))

    rep = check_field(th, a, b, c, p=0.5 + 5e-7, law="linear")
    dump("Williams with p=0.5+5e-7 (exponent tol 1e-6: should PASS exponent)", rep)
    results.append(("p off 5e-7", failed(rep)))

    # 7. Mode I field declared mode II: should fail symmetry
    rep = check_field(*williams_crack(), p=0.5, law="linear", mode="II")
    dump("Mode I Williams declared mode II", rep)
    results.append(("I as II", failed(rep)))

    # 8. Power-law field constructed from Airy with WRONG p but Airy-derived so
    # equilibrium ALWAYS holds for any p if stresses come from Airy.
    # Use Williams angular F (which satisfies faces for p=0.5 / lam=0.5) but
    # recompute stresses with a different p. Faces: Stt=(2-p)(1-p)F, Srt=-(1-p)F'
    # At faces F=F'=0 for the Williams F at +-pi, so faces still free for ANY p.
    # Equilibrium holds identically from Airy. Compatibility will fail unless p=0.5.
    # Exponent check will fail if p != 0.5.
    th = np.linspace(-np.pi, np.pi, 4001)
    lam = 0.5
    B = -np.cos((lam + 1) * np.pi) / np.cos((lam - 1) * np.pi)
    F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
    Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
    Fpp = -(lam + 1) ** 2 * np.cos((lam + 1) * th) - B * (lam - 1) ** 2 * np.cos(
        (lam - 1) * th
    )
    p_wrong = 0.3
    srr = Fpp + (2 - p_wrong) * F
    stt = (2 - p_wrong) * (1 - p_wrong) * F
    srt = -(1 - p_wrong) * Fp
    rep = check_field(th, srr, stt, srt, p=p_wrong, law="linear")
    dump("Airy-derived stresses with Williams F but p=0.3 (eq holds, not a solution)", rep)
    results.append(("Airy wrong p", failed(rep)))

    # 9. Same but declare p=0.5 with stresses built at p=0.3 — eq should fail
    rep = check_field(th, srr, stt, srt, p=0.5, law="linear")
    dump("Airy p=0.3 stresses declared p=0.5", rep)
    results.append(("mismatch p declare", failed(rep)))

    # 10. Sign error in shear (common polar-basis bug) — already in tests as 1.05
    th, a, b, c = williams_crack()
    rep = check_field(th, a, b, -c, p=0.5, law="linear")
    dump("Sign-flipped shear", rep)
    results.append(("flipped shear", failed(rep)))

    # 11. Smooth random-ish field that is traction-free and symmetric: F(theta)=(cos(th/2))**3 or similar
    # F even, F(pi)=0, F'(pi)=0. For crack: F(pi)=cos(pi/2)^3=0, F' ~ -1.5 cos^2 sin, at pi sin(pi/2)=1, cos=0 so F'=0.
    th = np.linspace(-np.pi, np.pi, 4001)
    F = np.cos(th / 2) ** 3
    Fp = -1.5 * np.cos(th / 2) ** 2 * np.sin(th / 2)
    Fpp = -1.5 * (
        -np.cos(th / 2) * np.sin(th / 2) ** 2 * 0.5 * 2
        + np.cos(th / 2) ** 2 * 0.5 * np.cos(th / 2)
    )
    # d/dth [cos^2 * sin] = 2 cos (-sin/2) sin + cos^2 (cos/2)
    Fpp = -1.5 * (
        2 * np.cos(th / 2) * (-0.5 * np.sin(th / 2)) * np.sin(th / 2)
        + np.cos(th / 2) ** 2 * (0.5 * np.cos(th / 2))
    )
    p = 0.5
    srr = Fpp + (2 - p) * F
    stt = (2 - p) * (1 - p) * F
    srt = -(1 - p) * Fp
    rep = check_field(th, srr, stt, srt, p=p, law="linear")
    dump("Fake even Airy F=cos^3(th/2), p=0.5 (eq+faces+symmetry+exponent hold; compatibility?)", rep)
    results.append(("fake Airy F", failed(rep)))

    # 12. Power-law: same fake Airy with HRR exponent n=3 p=0.25
    p = 0.25
    srr = Fpp + (2 - p) * F
    stt = (2 - p) * (1 - p) * F
    srt = -(1 - p) * Fp
    rep = check_field(th, srr, stt, srt, p=p, n=3, law="power", plane="strain")
    dump("Fake Airy F with HRR p=1/4 n=3 (eq, faces, energy, J-exponent hold)", rep)
    results.append(("fake Airy HRR p", failed(rep)))

    # 13. Linear combination of two eigenfields (not a pure eigenfunction)
    th, a, b, c = williams_crack()
    # add a little of T-stress
    srrT = np.cos(th) ** 2
    sttT = np.sin(th) ** 2
    srtT = -np.sin(th) * np.cos(th)
    rep = check_field(th, a + 0.2 * srrT, b + 0.2 * sttT, c + 0.2 * srtT, p=0.5, law="linear")
    dump("Williams + 0.2 T-stress, declared p=0.5 (not separable at one p)", rep)
    results.append(("Williams+T", failed(rep)))

    # 14. Plane strain vs plane stress: Williams is valid for both (linear).
    # Already checked.

    # 15. Notch with crack exponent p=0.5
    alpha = 2 * np.pi / 3
    thn, an, bn, cn, ptrue = williams_notch(alpha)
    rep = check_field(thn, an, bn, cn, p=0.5, law="linear", alpha=alpha)
    dump(f"True notch field (p={ptrue:.4f}) declared p=0.5 as if a crack", rep)
    results.append(("notch as crack p", failed(rep)))

    # 16. Crack field submitted as a notch (truncate theta, keep p=0.5)
    th, a, b, c = williams_crack()
    alpha = 2 * np.pi / 3
    keep = np.abs(th) <= alpha + 1e-12
    rep = check_field(th[keep], a[keep], b[keep], c[keep], p=0.5, law="linear", alpha=alpha)
    dump("Crack Williams truncated to 120deg notch, p=0.5", rep)
    results.append(("crack as notch", failed(rep)))

    # 17. Finite energy borderline: q-p+2 == 0  (p=2/(n+1) for power)
    # Need a field; use fake Airy. p=2/4=0.5, n=3. Energy: q=-1.5, q-p+2=0, passed = (0>0) False
    p = 2.0 / 4.0
    srr = Fpp + (2 - p) * F
    stt = (2 - p) * (1 - p) * F
    srt = -(1 - p) * Fp
    rep = check_field(th, srr, stt, srt, p=p, n=3, law="power", plane="strain")
    dump("Energy exactly zero exponent (p(n+1)=2): should FAIL finite energy", rep)
    results.append(("energy =0", failed(rep)))

    # 18. Reported J match on a notch: applies=False so won't fail even if reported_J is garbage
    thn, an, bn, cn, p = williams_notch(2 * np.pi / 3)
    rep = check_field(
        thn, an, bn, cn, p=p, law="linear", alpha=2 * np.pi / 3, reported_J=12345.0
    )
    dump("Notch with absurd reported_J=12345 (should not FAIL that check)", rep)
    r = by_name(rep).get("reported J (or I) matches the field's own J")
    print("  reported-J check:", r)
    results.append(("notch reported_J", failed(rep)))

    # 19. Half-range without enough points near face
    th = np.linspace(0, np.pi, 21)
    k = 1 / np.sqrt(2 * np.pi)
    a = k * (1.25 * np.cos(th / 2) - 0.25 * np.cos(1.5 * th))
    b = k * (0.75 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th))
    c = k * (0.25 * np.sin(th / 2) + 0.25 * np.sin(1.5 * th))
    rep = check_field(th, a, b, c, p=0.5, law="linear", mode="I")
    dump("Williams half-range only 21 points", rep)
    results.append(("half 21 pts", failed(rep)))

    # 20. Cartesian eq samples at r=1 only; a field that is equilibrated on the unit circle
    # but not as a separable field? Hard. Skip.

    # 21. nu ignored for power law: Williams power n=1 with plane stress
    th, a, b, c = williams_crack()
    rep = check_field(th, a, b, c, p=0.5, n=1, law="power", plane="stress")
    dump("Williams as n=1 power plane STRESS (szz=0 vs incompressible)", rep)
    results.append(("n=1 power PS", failed(rep)))

    # 22. Check Saint-Venant tolerance: 10*max(tol,1e-4)=1e-3 for default tol=1e-3
    print("\nSaint-Venant default tol =", max(1e-3, 1e-4) * 10)

    # 23. Linear notch Williams exponent check uses mode==I; mixed/II skipped
    # already.

    # 24. Wrong opening: use F for 120deg but alpha=135
    alpha_true = 2 * np.pi / 3
    alpha_wrong = 0.75 * np.pi
    thn, an, bn, cn, p = williams_notch(alpha_true)
    # re-sample? theta covers [0, alpha_true] not alpha_wrong — ValueError
    try:
        rep = check_field(thn, an, bn, cn, p=p, law="linear", alpha=alpha_wrong)
        dump("theta range mismatch", rep)
    except Exception as e:
        print("\n=== theta vs alpha mismatch ===")
        print(type(e).__name__, e)
        results.append(("theta mismatch", {type(e).__name__}))

    # 25. Compatibility: linear field, declare law=power n=1 plane strain — expect pass
    # already.

    print("\n\n######## SUMMARY ########")
    for item in results:
        print(item)


if __name__ == "__main__":
    main()
