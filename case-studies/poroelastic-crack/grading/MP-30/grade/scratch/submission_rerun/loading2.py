#!/usr/bin/env python3
"""
Loading 2: sigma_n = 0, p = H(t) on the faces.

K2: The *total* traction on the faces is zero, so there is no mode-I driving
force -> K2_0 = K2_inf = 0.

V2_inf: long-time (fully drained, steady pressure p=1 on the faces).
The crack opens because the pore pressure acts on the solid skeleton through
the Biot effective stress.  We derive V2_inf from the drained elastic crack
compliance acted on by the effective opening traction.

Derivation:
-----------
Effective stress (Biot):  sigma'_ij = sigma_ij + alpha p delta_ij.
On the crack faces (exterior total traction sigma_n = 0, pressure p = 1):
The net traction that opens the crack is the *effective* face traction
    sigma'_n = sigma_n + alpha p = 0 + alpha = alpha  (tension positive,
    i.e. the skeleton is pulled open by alpha*p).

So the opening is that of an elastic crack (drained properties nu, G) under a
uniform effective tension alpha:

    V2_inf = 2 pi (1 - nu)/G * (alpha * 1) * L^2

This is the standard poroelastic result: a crack inflated by pore pressure
p opens with an equivalent elastic pressure alpha p (for the drained solid).

Small time: at t=0+ the pressure is still p=1 on the faces but nothing has
diffused; the instantaneous (undrained) elastic response would give
2 pi (1-nu_u)/G alpha L^2.  The problem only asks for V2_inf.
"""
import numpy as np


def loading2(nu, nuu, alpha, G, L=1.0):
    K0 = 0.0
    Kinf = 0.0
    Vinf = 2.0 * np.pi * (1.0 - nu) / G * alpha * L * L
    return K0, Kinf, Vinf


if __name__ == "__main__":
    nu, nuu, alpha, G, L = 0.2, 0.4, 0.89, 1.0, 1.0
    print(loading2(nu, nuu, alpha, G, L))