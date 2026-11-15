import numpy as np
from fem import solve_fem

print("DCB check at fine settings (psi must be ~0):")
for dmin in [4e-5, 1.5e-5, 6e-6]:
    r = solve_fem(1.0, 1.0, 60.0, 1.0, 1.0, 1/3, 0.0, 1.0, dmin=dmin,
                  ny1=18, ny2=18, fit_rmin=8 * dmin, fit_rmax=0.15, nterms=3)
    print(f"  dmin={dmin:g}: |K|={r['Kabs']:.5f} (err {abs(r['Kabs']/np.sqrt(12)-1)*100:+.2f}%) "
          f"psi={r['psi_deg']:+.4f} ndof={r['ndof']}")

print("target B, ultra-fine:")
for dmin in [4e-5, 1.5e-5, 6e-6]:
    r = solve_fem(1.0, 2.0, 60.0, 9.0, 1.0, 1/3, 0.0, 1.0, dmin=dmin,
                  ny1=18, ny2=28, fit_rmin=8 * dmin, fit_rmax=0.15, nterms=3)
    print(f"  dmin={dmin:g}: |K|={r['Kabs']:.5f} psi={r['psi_deg']:+.4f} ndof={r['ndof']}")
print("target A, ultra-fine:")
for dmin in [4e-5, 1.5e-5, 6e-6]:
    r = solve_fem(1.0, 2.0, 60.0, 9.0, 1.0, 1/3, 1.0, 0.0, dmin=dmin,
                  ny1=18, ny2=28, fit_rmin=8 * dmin, fit_rmax=0.15, nterms=3)
    print(f"  dmin={dmin:g}: |K|={r['Kabs']:.5f} psi={r['psi_deg']:+.4f} ndof={r['ndof']}")
print("target AB, ultra-fine:")
for dmin in [4e-5, 1.5e-5, 6e-6]:
    r = solve_fem(1.0, 2.0, 60.0, 9.0, 1.0, 1/3, 1.0, 1.0, dmin=dmin,
                  ny1=18, ny2=28, fit_rmin=8 * dmin, fit_rmax=0.15, nterms=3)
    print(f"  dmin={dmin:g}: |K|={r['Kabs']:.5f} psi={r['psi_deg']:+.4f} ndof={r['ndof']}")
