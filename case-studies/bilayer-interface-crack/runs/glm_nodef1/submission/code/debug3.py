import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import fit_series

nu = 1/3
D = [plane_strain_D(1.0, nu)]*2
h = H = 1.0
L = 30.0
P, M = 0.0, 1.0
sigma = {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
         2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}

for (d1, d2) in [(5e-5, 2e-4), (2.5e-5, 1e-4)]:
    mesh = BilayerCrack(h, H, L, d1, d2, q1=1.35, q2=1.5, cap2=0.35)
    f = end_loads(mesh, sigma)
    u = solve(mesh, D, f)
    rs = np.geomspace(5e-4, 0.2, 40)
    T = interface_tractions(mesh, D, u, rs)
    y = T*np.sqrt(2*np.pi*rs)
    print(f"mesh d1={d1:.0e} d2={d2:.0e} ndof={2*len(mesh.coords)}")
    for rmax in [0.01, 0.02, 0.05]:
        for rmin in [5e-4, 1e-3, 3e-3]:
            m = (rs >= rmin) & (rs <= rmax)
            for pw in [(0.0, 0.5, 1.0), (0.0, 0.5, 1.0, 1.5)]:
                c = fit_series(rs[m], y[m], pw)
                print(f"   rmin={rmin:.0e} rmax={rmax:.3f} pw={len(pw)}: K={c[0].real:8.4f}"
                      + (f"   <<< {len(pw)} terms" if len(pw) == 3 else ""))
    rd, dd = crack_opening(mesh, u)
    # delta fit: delta = i*C*K*sqrt(r)*(1+c sqrt r), C = 4*(2/Ebar) = 64/9
    C = 4*(2*8/9)
    m = (rd > 5e-4) & (rd < 0.2)
    yd = dd[m]/np.sqrt(rd[m])/ (1j*C)
    for rmax in [0.01, 0.02, 0.05]:
        mm = (rd[m] <= rmax)
        c = fit_series(rd[m][mm], yd[mm], (0.0, 0.5, 1.0))
        print(f"   delta-fit rmin=5e-4 rmax={rmax:.3f}: K={c[0].real:8.4f}")
print("exact K_I = 3.464102")
