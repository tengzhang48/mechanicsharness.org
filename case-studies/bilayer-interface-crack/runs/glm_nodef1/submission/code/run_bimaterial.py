"""Bimaterial runs v2: cases A and B with fixed complex fits + diagnostics."""
import numpy as np
import sys, json
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import (fit_series, measure_eps, measure_eps_delta,
                  K_from_traction, amplitude_from_delta, measure_eps_unwrap)

nu = 1/3
E1, E2 = 9.0, 1.0
h, H = 1.0, 2.0
L = 30.0
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]
Eb1, Eb2 = E1/(1-nu**2), E2/(1-nu**2)
eps_an = np.log((1-0.2)/(1+0.2))/(2*np.pi)     # = ln(2/3)/(2pi) = -0.0645033
G_A = 1/(2*Eb1*h) + 1/(2*Eb2*H) + 6*((h+H)/2)**2/(Eb2*H**3)
G_B = 6*(1/(Eb1*h**3) + 1/(Eb2*H**3))
S = 1/Eb1 + 1/Eb2
print(f"eps_an = {eps_an:+.6f};  G_A E2 h/P^2 = {G_A:.7f};  G_B E2 h^3/M^2 = {G_B:.7f}")
print(f"G-K candidates E2*S/2*{{1, sech^2, cosh^2}} = "
      f"{E2*S/2:.6f}, {E2*S/2/np.cosh(np.pi*eps_an)**2:.6f}, {E2*S/2*np.cosh(np.pi*eps_an)**2:.6f}")


def sigma(P, M):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}


def run_case(P, M, d1, d2, rmin=3e-3, rmax=6e-2, Lx=L):
    mesh = BilayerCrack(h, H, Lx, d1, d2, q1=1.35, q2=1.5, cap2=0.35)
    f = end_loads(mesh, sigma(P, M))
    u = solve(mesh, D, f)
    rs = np.geomspace(rmin, rmax, 40)
    T = interface_tractions(mesh, D, u, rs)
    rd, dd = crack_opening(mesh, u)
    m = (rd >= rmin) & (rd <= rmax)
    return dict(rs=rs, T=T, rd=rd[m], dd=dd[m], ndof=2*len(mesh.coords))


def analyse(R, eps):
    K = K_from_traction(R["rs"], R["T"], eps)
    Ad = amplitude_from_delta(R["rd"], R["dd"], -eps)
    return K, Ad


out = {}
for case, (P, M) in [("A", (1.0, 0.0)), ("B", (0.0, 1.0))]:
    print(f"\n========== case {case}")
    rows = []
    for (name, d1, d2) in [("L1", 6e-5, 1.5e-4), ("L2", 3e-5, 7.5e-5),
                           ("L3", 1.5e-5, 3.75e-5)]:
        R = run_case(P, M, d1, d2)
        eps_scan, resT, _ = measure_eps(R["rs"], R["T"])
        eps_unw, _ = measure_eps_unwrap(R["rs"], R["T"]*np.sqrt(2*np.pi*R["rs"]))
        epsd_scan, resd, _ = measure_eps_delta(R["rd"], R["dd"])
        # delta slope measured directly too
        epsd_unw, _ = measure_eps_unwrap(R["rd"], R["dd"]/np.sqrt(R["rd"]))
        Kp, Ad = analyse(R, eps_unw)
        rows.append(dict(name=name, ndof=R["ndof"], eps_scan=eps_scan, eps_unw=eps_unw,
                         resT=resT, epsd_scan=epsd_scan, epsd_unw=epsd_unw, resd=resd,
                         K=Kp, Ad=Ad))
        print(f" {name} ndof={R['ndof']}: eps_T scan={eps_scan:+.5f} unw={eps_unw:+.5f} "
              f"(res {resT:.1e}); eps_d scan={epsd_scan:+.5f} unw={epsd_unw:+.5f}")
        print(f"      K = {Kp:.6f}   |K|={abs(Kp):.6f}  phase={np.degrees(np.angle(Kp)):+.4f} deg")
        print(f"      A_delta = {Ad:.6f}")
    out[case] = rows
    K = rows[-1]["K"]
    G = G_A if case == "A" else G_B
    print(f"  E2*G/|K|^2 = {E2*G/abs(K)**2:.6f}  (sech^2 candidate {E2*S/2/np.cosh(np.pi*eps_an)**2:.6f})")

json.dump({c: [{k: complex(v) if isinstance(v, (complex, np.complex128)) else float(v)
                for k, v in r.items()} for r in v] for c, v in out.items()},
          open('/work/code/bimaterial_raw2.json', 'w'), indent=1)
print("saved")
