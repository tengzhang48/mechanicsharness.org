#!/usr/bin/env python3
"""
Exact mode-I (uniform-pressure crack) stress field on y=0, to fix sigma_xx at
the crack face.  This determines the undrained pore pressure next to the face
for loading 1 (t=0+).

Westergaard: for a crack -L<=x<=L, uniform pressure P (opening), the stress
field on the continuation of the crack line and on the faces:

  zeta = z + sqrt(z^2 - L^2),  Phi(z) = P (z/sqrt(z^2-L^2) - 1) [for opening
  traction], classic:

  sigma_xx = Re Phi - y Im Phi'
  sigma_yy = Re Phi + y Im Phi'
  sigma_xy = -y Re Phi'   (signs per convention)

The traction on the crack faces (|x|<L, y=0) from the material is
  sigma_yy = -P  (so applied opening pressure P balances it).
And on the faces, the in-plane stress sigma_xx = ?
"""
import numpy as np


def sigma_on_crackline(x, L, P, y=1e-9):
    """Stress components at point (x, y) using Westergaard for a uniformly
    pressurized crack (opening traction P on faces)."""
    z = x + 1j * y
    # Phi(z) for uniform pressure: Phi(z) = P * ( z/sqrt(z^2-L^2) - 1 )
    root = np.sqrt(z * z - L * L)
    Phi = P * (z / root - 1.0)
    dPhi = P * (-L * L) / root**3  # d/dz [z/root] = -L^2/root^3
    sxx = Phi.real - y * dPhi.imag
    syy = Phi.real + y * dPhi.imag
    sxy = -y * dPhi.real
    return sxx, syy, sxy


if __name__ == "__main__":
    L, P = 1.0, 1.0
    print("On crack faces (|x|<L):  expect sigma_yy -> -P, sigma_xy -> 0")
    for x in (0.0, 0.3, 0.7, 0.95):
        for y in (1e-6,):
            sxx, syy, sxy = sigma_on_crackline(x, L, P, y)
            print(f"  x={x:4.2f} y={y:.0e}: sxx={sxx: .6f} syy={syy: .6f} sxy={sxy: .2e}")
    print("Just outside tip (x>L): means the K-dominated field")
    for x in (1.001, 1.01, 1.1):
        sxx, syy, sxy = sigma_on_crackline(x, L, P, 0.0)
        print(f"  x={x:5.3f}: sxx={sxx: .6f} syy={syy: .6f}")
    print("Toward the face interior limit (x<L, y->0+):")
    for x in (0.5,):
        for y in (1e-3, 1e-6, 1e-9):
            sxx, syy, sxy = sigma_on_crackline(x, L, P, y)
            print(f"  x={x} y={y:.0e}: sxx={sxx: .6f}")