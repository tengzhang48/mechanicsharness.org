"""Validation 2: finite interface crack |x1| < a on y = 0 between two bonded half-planes (material 1 above),
remote sigma22 = sigma, sigma12 = tau (and the compatible sigma11 in material 2), in a square [-W,W]^2 whose
boundary carries the uncracked remote tractions.  Exact right-tip SIF (Rice & Sih 1965 / Rice 1988):
    K = (sigma + i tau)(1 + 2 i eps) sqrt(pi a) (2a)^(-i eps)."""
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
from bilayer import RINGS


def run(E1, nu1, E2, nu2, W, m=1.0, a=1.0, R0=0.4, rmin=1e-6, tag=""):
    nh = int(round(4 * m))
    pts, cells, mat, info = build_mesh(-W, W, -W, W, [(-a, 1), (a, -1)], [(-a, a)], R0=R0, nh=nh, beta=0.15 / m,
                                       smax_x=2.0 / m, smax_y=0.1 * W / m, rmin=rmin)
    chk = check_mesh(pts, cells, -W, W, -W, W, 2 * a)
    assert chk["min_area"] > 0 and chk["area_err"] < 1e-12 and chk["blen_err"] < 1e-12 and chk["unused_nodes"] == 0
    mu1, mu2 = E1 / (2 * (1 + nu1)), E2 / (2 * (1 + nu2))

    def loads(sig, tau):
        s11_1 = 0.0
        s11_2 = sig * (nu2 - nu1 * mu2 / mu1) / (1 - nu2)
        tol = 1e-9 * W
        return [
            (lambda xm: np.abs(xm[:, 1] - W) < tol, lambda x: ufl.as_vector((tau, sig))),
            (lambda xm: np.abs(xm[:, 1] + W) < tol, lambda x: ufl.as_vector((-tau, -sig))),
            (lambda xm: (np.abs(xm[:, 0] - W) < tol) & (xm[:, 1] > 0), lambda x: ufl.as_vector((s11_1, tau))),
            (lambda xm: (np.abs(xm[:, 0] - W) < tol) & (xm[:, 1] < 0), lambda x: ufl.as_vector((s11_2, tau))),
            (lambda xm: (np.abs(xm[:, 0] + W) < tol) & (xm[:, 1] > 0), lambda x: ufl.as_vector((-s11_1, -tau))),
            (lambda xm: (np.abs(xm[:, 0] + W) < tol) & (xm[:, 1] < 0), lambda x: ufl.as_vector((-s11_2, -tau))),
        ]

    cases = {"tension": (1.0, 0.0), "shear": (0.0, 1.0)}
    pins = [((-W, -W), 0), ((-W, -W), 1), ((W, -W), 1)]
    t0 = time.time()
    sols = fe.solve(pts, cells, {1: (E1, nu1), 2: (E2, nu2)}, [loads(*c) for c in cases.values()], pins)
    res = dict(W=W, m=m, ndofs=sols[0].ndofs, time=time.time() - t0, cases={})
    for (cname, (sig, tau)), S in zip(cases.items(), sols):
        ring = fe.interaction_integrals(S, (a, 0.0), RINGS)
        eps = S.eps
        Kex = (sig + 1j * tau) * (1 + 2j * eps) * np.sqrt(np.pi * a) * (2 * a) ** (-1j * eps)
        Kalt = (sig + 1j * tau) * (1 - 2j * eps) * np.sqrt(np.pi * a) * (2 * a) ** (1j * eps)  # opposite-eps convention
        Gex = S.Gfac * abs(Kex) ** 2
        rc, tc = fe.consistent_tractions(S, a, a + 0.5)
        r, tu, tl = fe.interface_tractions(S, (a, 0.0), 0.5)
        Kc, _ = fe.traction_fit(rc, tc, eps, 3e-4, 3e-2, 3)
        Kn, _ = fe.traction_fit(r, 0.5 * (tu + tl), eps, 3e-4, 3e-2, 3)
        out = dict(eps=eps, K_exact=[Kex.real, Kex.imag], psi_exact=float(np.degrees(np.angle(Kex))),
                   psi_opposite_convention=float(np.degrees(np.angle(Kalt))),
                   rings=[dict(r1=d["r1"], r2=d["r2"], K=[d["K"].real, d["K"].imag],
                               absK_relerr=abs(d["K"]) / abs(Kex) - 1,
                               dpsi_deg=float(np.degrees(np.angle(d["K"] / Kex))), J_over_Gex=d["J"] / Gex)
                          for d in ring],
                   consistent_traction_fit=dict(K=[Kc.real, Kc.imag], absK_relerr=abs(Kc) / abs(Kex) - 1,
                                                dpsi_deg=float(np.degrees(np.angle(Kc / Kex)))),
                   nodal_traction_fit=dict(K=[Kn.real, Kn.imag], absK_relerr=abs(Kn) / abs(Kex) - 1,
                                           dpsi_deg=float(np.degrees(np.angle(Kn / Kex)))))
        res["cases"][cname] = out
        print(f"[{tag} W={W} m={m} ndofs={res['ndofs']} t={res['time']:.1f}s] {cname}: eps={eps:.6f} "
              f"psi_exact={out['psi_exact']:.4f} (opposite-eps convention would give {out['psi_opposite_convention']:.4f})")
        for d in out["rings"]:
            print(f"    M-int ring {d['r1']:.2f}-{d['r2']:.2f}: |K| rel err={d['absK_relerr']:+.2e}  "
                  f"psi err={d['dpsi_deg']:+.4f} deg  J/G_exact-1={d['J_over_Gex'] - 1:+.2e}")
        for nm in ("consistent_traction_fit", "nodal_traction_fit"):
            d = out[nm]
            print(f"    {nm}: |K| rel err={d['absK_relerr']:+.2e}  psi err={d['dpsi_deg']:+.4f} deg")
    return res


if __name__ == "__main__":
    W = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    m = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    which = sys.argv[3] if len(sys.argv) > 3 else "S2"
    mats = {"S2": (9.0, 1 / 3, 1.0, 1 / 3), "S1": (5.0, 0.22, 1.0, 0.35), "S1swap": (1.0, 0.35, 5.0, 0.22)}[which]
    res = run(*mats, W=W, m=m, tag=which)
    os.makedirs("results", exist_ok=True)
    with open(f"results/griffith_{which}_W{W:g}_m{m:g}.json", "w") as f:
        json.dump(res, f, indent=1)
