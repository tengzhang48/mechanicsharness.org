"""Check the energy relation G = (c1+c2)|K|^2/(16 cosh^2(pi eps)) on the auxiliary field itself:
contour J-integral of the pure asymptotic field around a circle, for several Ka (and two radii),
and the crack-opening relation  d2 + i d1 = 8 K r^(i eps) sqrt(r/2pi) / ((1+2i eps) cosh(pi eps) E*)."""
import numpy as np
from asym import TipField
from numpy.polynomial.legendre import leggauss

for (E1, nu1, E2, nu2) in [(5.0, 0.22, 1.0, 0.35), (9.0, 1/3, 1.0, 1/3), (1.0, 0.3, 1.0, 0.3)]:
    tf = TipField(E1, nu1, E2, nu2)
    c1 = (tf.k1 + 1) / tf.mu1; c2 = (tf.k2 + 1) / tf.mu2
    fac = (c1 + c2) / (16 * np.cosh(np.pi * tf.eps) ** 2)
    xg, wg = leggauss(200)
    for Ka in [1.0, 1j, 0.6 - 0.8j]:
        for R in [0.3, 2.0]:
            J = 0.0
            for (t0, t1) in [(-np.pi, 0.0), (0.0, np.pi)]:
                t = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0); w = 0.5 * (t1 - t0) * wg * R
                x, y = R * np.cos(t), R * np.sin(t)
                y = np.where(t >= 0, np.maximum(y, 1e-300), np.minimum(y, -1e-300))
                F = tf.field(x, y, Ka)
                n1, n2 = np.cos(t), np.sin(t)
                W = 0.5 * (F["s11"] * F["u11"] + F["s22"] * F["u22"] + F["s12"] * (F["u12"] + F["u21"]))
                J += np.sum(w * (W * n1 - (F["s11"] * n1 + F["s12"] * n2) * F["u11"] - (F["s12"] * n1 + F["s22"] * n2) * F["u21"]))
            print(f"E1/E2={E1/E2:g} Ka={Ka} R={R}: J/(fac|Ka|^2) - 1 = {J/(fac*abs(Ka)**2)-1:.2e}")
    # crack-opening displacement relation
    Ebar1 = E1/(1-nu1**2); Ebar2 = E2/(1-nu2**2); Estar = 2/(1/Ebar1 + 1/Ebar2)
    Ka = 0.6 - 0.8j; r = np.array([0.01, 0.1, 1.0])
    Fu = tf.field(-r, 0*r, Ka); Fl = tf.field(-r, -1e-300*np.ones_like(r), Ka)
    d = (Fu["u2"] - Fl["u2"]) + 1j*(Fu["u1"] - Fl["u1"])
    ref = 8*Ka*r**(1j*tf.eps)*np.sqrt(r/(2*np.pi))/((1+2j*tf.eps)*np.cosh(np.pi*tf.eps)*Estar)
    print("   COD relation max rel err:", np.max(np.abs(d/ref - 1)))
