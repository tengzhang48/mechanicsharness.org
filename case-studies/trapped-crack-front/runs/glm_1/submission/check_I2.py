"""Validate corrected I1/I2 for b = e^{-X}:
  PV int e^{-X}/(Xj-X) dX = J = e^{-Xj} Ei(Xj)
  fp  int e^{-X}/(Xj-X)^2 dX = J - 1/Xj
"""
import numpy as np
from scipy.special import expi
from kernels_x import BasisFun, Jfun, nodes_for

bf = BasisFun(0.0, 0)
for Xj in [0.2, 0.5, 1.0, 3.0]:
    X, W = nodes_for(Xj)
    B = bf.eval(X)
    b0, d1, d2 = bf.derivs(Xj)
    t = Xj - X; u = X - Xj
    phi = np.exp(-u)
    I1 = ((B - b0*phi)/t) @ W + b0*expi(Xj)
    ph0 = np.exp(-u)*(1 + u + 0.5*u*u)
    ph1 = np.exp(-u)*(u + u*u)
    ph2 = np.exp(-u)*0.5*u*u
    S = b0*ph0 + d1*ph1 + d2*ph2
    ex = np.exp(Xj); Ei = expi(Xj)
    I2 = ((B - S)/t**2) @ W + b0*ex*(0.5 - 1.0/Xj) + d1*(-Ei + ex) + d2*0.5*ex
    Jex = Jfun(Xj)
    print(f"Xj={Xj}: I1={I1:.10f} exact={Jex:.10f} | I2={I2:.10f} exact={Jex - 1/Xj:.10f}")
