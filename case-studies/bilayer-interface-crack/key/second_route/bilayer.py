"""Bilayer interface crack (semi-infinite crack on x<0, y=0), loads applied at x = -Lb as the beam stresses.
Units: h = 1, E2 = 1, P = 1 (case A) or M = 1 (case B)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
import json
import time
import numpy as np
import ufl
from meshgen import build_mesh, check_mesh
import fe

SETS = {
    "V1": (1.0, 0.3, 1.0, 0.3, 1.0),      # validation: homogeneous, h = H
    "S1": (5.0, 0.22, 1.0, 0.35, 0.3),
    "S2": (9.0, 1 / 3, 1.0, 1 / 3, 0.5),
    "S3": (1.0, 0.3, 1.0, 0.3, 0.5),
    "S4": (1.0, 0.3, 1.0, 0.3, 0.1),
}

RINGS = [(0.05, 0.25), (0.1, 0.4), (0.2, 0.6), (0.3, 0.85)]


def beam_G(E1, nu1, E2, nu2, h, H, P, M):
    Eb1, Eb2 = E1 / (1 - nu1 ** 2), E2 / (1 - nu2 ** 2)
    return (P ** 2 / h + 12 * M ** 2 / h ** 3) / (2 * Eb1) + (P ** 2 / H + 12 * (M + P * (h + H) / 2) ** 2 / H ** 3) / (2 * Eb2)


def run(setname, m=1.0, Lfac=4.0, R0=0.4, rmin=1e-6, verbose=True):
    E1, nu1, E2, nu2, hH = SETS[setname]
    h, H = 1.0, 1.0 / hH
    L = Lfac * (h + H)
    X0, X1 = -L, L
    nh = int(round(4 * m))
    beta = 0.15 / m
    smax_x = min(1.0, 0.15 * (h + H)) / m
    smax_y = 0.15 * (h + H) / m
    t0 = time.time()
    pts, cells, mat, info = build_mesh(X0, X1, -H, h, [(0.0, -1)], [(X0 - 1.0, 0.0)], R0=R0, nh=nh, beta=beta,
                                       smax_x=smax_x, smax_y=smax_y, rmin=rmin)
    chk = check_mesh(pts, cells, X0, X1, -H, h, L)
    assert chk["min_area"] > 0 and chk["area_err"] < 1e-12 and chk["blen_err"] < 1e-12 and chk["unused_nodes"] == 0

    def loads(P, M):
        left_up = lambda xm: (np.abs(xm[:, 0] - X0) < 1e-9) & (xm[:, 1] > 0)
        left_lo = lambda xm: (np.abs(xm[:, 0] - X0) < 1e-9) & (xm[:, 1] < 0)
        s1 = lambda x: -P / h - 12 * M * (x[1] - h / 2) / h ** 3
        s2 = lambda x: P / H + 12 * (M + P * (h + H) / 2) * (x[1] + H / 2) / H ** 3
        # outward normal (-1, 0): traction = (-s11, 0)
        return [(left_up, lambda x: ufl.as_vector((-s1(x), 0.0))),
                (left_lo, lambda x: ufl.as_vector((-s2(x), 0.0)))]

    cases = {"A": (1.0, 0.0), "B": (0.0, 1.0)}
    pins = [((X1, -H), 0), ((X1, -H), 1), ((X1, h), 0)]
    sols = fe.solve(pts, cells, {1: (E1, nu1), 2: (E2, nu2)}, [loads(*cases["A"]), loads(*cases["B"])], pins)
    t1 = time.time()
    res = dict(set=setname, m=m, Lfac=Lfac, R0=R0, ndofs=sols[0].ndofs, ncells=len(cells), Nr=info["Nr"],
               Nt=info["Nt"], max_aspect=float(chk["max_aspect"]), solve_time=t1 - t0, cases={})
    for (cname, (P, M)), S in zip(cases.items(), sols):
        ring = fe.interaction_integrals(S, (0.0, 0.0), RINGS)
        eps = S.eps
        G = beam_G(E1, nu1, E2, nu2, h, H, P, M)
        KG = np.sqrt(G / S.Gfac)
        r, tu, tl = fe.interface_tractions(S, (0.0, 0.0), 0.5)
        tav = 0.5 * (tu + tl)
        fits = {}
        for (rlo, rhi) in [(1e-4, 1e-2), (3e-4, 3e-2), (1e-3, 0.1)]:
            Kf, n = fe.traction_fit(r, tav, eps, rlo, rhi, 3)
            fits[f"{rlo:g}-{rhi:g}"] = dict(K=[Kf.real, Kf.imag], absK=abs(Kf),
                                            psi=float(np.degrees(np.angle(Kf * h ** (1j * eps)))), n=n)
        rc, tc = fe.consistent_tractions(S, 0.0, 0.5)
        cfits = {}
        for (rlo, rhi) in [(1e-4, 1e-2), (3e-4, 3e-2), (1e-3, 0.1)]:
            Kf, n = fe.traction_fit(rc, tc, eps, rlo, rhi, 3)
            cfits[f"{rlo:g}-{rhi:g}"] = dict(K=[Kf.real, Kf.imag], absK=abs(Kf),
                                             psi=float(np.degrees(np.angle(Kf * h ** (1j * eps)))), n=n)
        out = dict(eps=eps, G_beam=G, absK_from_G=KG, resid=S.resid, consistent_traction_fit=cfits,
                   rings=[dict(r1=d["r1"], r2=d["r2"], K=[d["K"].real, d["K"].imag], absK=abs(d["K"]),
                               psi=float(np.degrees(np.angle(d["K"] * h ** (1j * eps)))), J=d["J"],
                               J_over_G=d["J"] / G) for d in ring],
                   traction_fit=fits)
        res["cases"][cname] = out
        if verbose:
            print(f"[{setname} m={m} Lfac={Lfac} ndofs={res['ndofs']} t={t1 - t0:.1f}s] case {cname}: eps={eps:.6f} "
                  f"|K|_G={KG:.6f} resid={S.resid:.1e}")
            for d in out["rings"]:
                print(f"    M-int ring {d['r1']:.2f}-{d['r2']:.2f}: psi={d['psi']:.4f}  |K|={d['absK']:.6f} "
                      f"(|K|/|K|_G-1={d['absK'] / KG - 1:+.2e})  J/G-1={d['J_over_G'] - 1:+.2e}")
            for k, d in fits.items():
                print(f"    traction fit r in {k}: psi={d['psi']:.4f} |K|={d['absK']:.6f} "
                      f"(|K|/|K|_G-1={d['absK'] / KG - 1:+.2e}) n={d['n']}")
            for k, d in cfits.items():
                print(f"    consistent-traction fit r in {k}: psi={d['psi']:.4f} |K|={d['absK']:.6f} "
                      f"(|K|/|K|_G-1={d['absK'] / KG - 1:+.2e}) n={d['n']}")
    return res


if __name__ == "__main__":
    setname = sys.argv[1]
    m = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    Lfac = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
    res = run(setname, m, Lfac)
    os.makedirs("results", exist_ok=True)
    with open(f"results/{setname}_m{m:g}_L{Lfac:g}.json", "w") as f:
        json.dump(res, f, indent=1)
