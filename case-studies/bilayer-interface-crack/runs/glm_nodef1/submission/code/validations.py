"""Final validations: homogeneous case A, combined load, L sensitivity, extra level."""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import fit_series

nu = 1/3


def make_D(E1, E2):
    return [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]


def sigma_of(P, M, h, H):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}


def Kfit(P, M, h, H, D, d1, d2, L=30.0, eps=0.0, windows=((5e-4, 1e-2), (1e-3, 2e-2), (2e-3, 4e-2))):
    mesh = BilayerCrack(h, H, L, d1, d2, q1=1.35, q2=1.5, cap2=0.35)
    f = end_loads(mesh, sigma_of(P, M, h, H), )
    u = solve(mesh, D, f)
    rs = np.geomspace(4e-4, 8e-2, 60)
    T = interface_tractions(mesh, D, u, rs)
    out = []
    for (rmin, rmax) in windows:
        m = (rs >= rmin) & (rs <= rmax)
        y = T[m]*np.sqrt(2*np.pi*rs[m])*np.exp(-1j*eps*np.log(rs[m]))
        c = fit_series(rs[m], y, (0.0, 0.5, 1.0))
        out.append(c[0])
    rd, dd = crack_opening(mesh, u)
    md = (rd >= 1e-3) & (rd <= 2e-2)
    yd = dd[md]/np.sqrt(rd[md])*np.exp(+1j*eps*np.log(rd[md]))
    cd = fit_series(rd[md], yd, (0.0, 0.5, 1.0))
    return out, cd[0], 2*len(mesh.coords)


print("=== homogeneous case A (P=1): exact |K| = 1, phase = +-90 deg")
D = make_D(1.0, 1.0)
for (d1, d2) in [(3e-5, 7.5e-5), (7.5e-6, 1.875e-5), (4e-6, 1e-5)]:
    out, Ad, ndof = Kfit(1.0, 0.0, 1.0, 1.0, D, d1, d2)
    for K in out:
        print(f"  d={d1:.1e}: K={K:.6f} |K|={abs(K):.6f} ph={np.degrees(np.angle(K)):+.3f}")
print("  (delta fit not usable for pure sliding sign check here)")

print("\n=== homogeneous case B (M=1): exact K = 2 sqrt(3), phase 0")
for (d1, d2) in [(7.5e-6, 1.875e-5)]:
    out, Ad, ndof = Kfit(0.0, 1.0, 1.0, 1.0, D, d1, d2)
    for K in out:
        print(f"  d={d1:.1e}: K={K:.6f} |K|={abs(K):.6f} ph={np.degrees(np.angle(K)):+.4f}"
              f"   err={abs(K)/(2*np.sqrt(3))-1:+.4%}")
    C = 4*(2*8/9)/np.sqrt(2*np.pi)
    print(f"  delta-fit K = {Ad/(1j*C):.6f}")

print("\n=== bimaterial: combined (P=1,M=1) direct vs superposition")
Db = make_D(9.0, 1.0)
eps = np.log(2/3)/(2*np.pi)
KA_list, KB_list = [], []
for (d1, d2) in [(3e-5, 7.5e-5), (7.5e-6, 1.875e-5), (4e-6, 1e-5)]:
    oA, dA, _ = Kfit(1.0, 0.0, 1.0, 2.0, Db, d1, d2, eps=eps)
    oB, dB, _ = Kfit(0.0, 1.0, 1.0, 2.0, Db, d1, d2, eps=eps)
    oC, dC, _ = Kfit(1.0, 1.0, 1.0, 2.0, Db, d1, d2, eps=eps)
    KA_list.append(oA[1]); KB_list.append(oB[1])
    print(f"  d={d1:.1e}: KA+KB = {oA[1]+oB[1]:.6f} (ph {np.degrees(np.angle(oA[1]+oB[1])):+.4f})  "
          f"direct = {oC[1]:.6f} (ph {np.degrees(np.angle(oC[1])):+.4f})")
print("  exact-relation magnitudes: |KA| = sqrt(GA/c):")
Eb1, Eb2 = 9/(1-nu**2), 1/(1-nu**2)
S = 1/Eb1+1/Eb2
cc = S/2/np.cosh(np.pi*eps)**2
GA = 1/(2*Eb1)+1/(2*Eb2*2)+6*(1.5)**2/(Eb2*8)
GB = 6*(1/Eb1+1/(Eb2*8))
print(f"   c={cc:.8f} |KA|={np.sqrt(GA/cc):.6f} |KB|={np.sqrt(GB/cc):.6f}")
print(f"   phase(KA+KB with exact mags) = "
      f"{np.degrees(np.angle(np.sqrt(GA/cc)*KA_list[-1]/abs(KA_list[-1]) + np.sqrt(GB/cc)*KB_list[-1]/abs(KB_list[-1]))):+.4f}")

print("\n=== L sensitivity (case B, d=7.5e-6):")
for L in [15.0, 30.0, 45.0]:
    out, Ad, ndof = Kfit(0.0, 1.0, 1.0, 2.0, Db, 7.5e-6, 1.875e-5, L=L, eps=eps)
    print(f"  L={L}: K={out[1]:.6f}")
