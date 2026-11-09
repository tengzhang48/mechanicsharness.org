import numpy as np
from scipy.optimize import brentq

E, nu, sy, J, Czz0 = 70000.0, 0.30, 10.0, 0.800, -0.005
A_hyd = (1 - 2 * nu) / (3 * E)
A_dev = lambda T: 1.5 / T + (nu - 0.5) / E
def props(T):
    G = 1 / (2 * A_dev(T)); K = 1 / (9 * A_hyd)
    Ee = 9 * K * G / (3 * K + G); ne = (3 * K - 2 * G) / (2 * (3 * K + G))
    return G, K, Ee, ne

# closed form: E_eff=T, nu_eff = 0.5 - 0.2 T/E, and nu_eff/sqrt(T) = -Czz0/2*sqrt(2pi/J)
alpha = -Czz0 / 2 * np.sqrt(2 * np.pi / J)
b = 1 / (5 * E)          # nu_eff = 0.5 - T/(5E)
x = (-alpha + np.sqrt(alpha**2 + 4 * b * 0.5)) / (2 * b)   # x = sqrt(T)
T_closed = x**2
T = brentq(lambda T: -2 * props(T)[3] / props(T)[2] * np.sqrt(J * props(T)[2] / (2 * np.pi)) - Czz0,
           1000, 60000, xtol=1e-13, rtol=1e-15)
G, K, Ee, ne = props(T)
KA = np.sqrt(J * Ee / (2 * np.pi))
print("T closed form = %.12f" % T_closed)
print("T numeric     = %.12f" % T)
print("E_eff = %.12f  nu_eff = %.12f  G_eff = %.12f  K_A = %.12f" % (Ee, ne, G, KA))

c30 = np.cos(np.pi/6)
Srr = KA * (1.25*c30 - 0.25*np.cos(np.pi/2))
Stt = KA * (0.75*c30 + 0.25*np.cos(np.pi/2))
Srt = KA * (0.25*np.sin(np.pi/6) + 0.25*np.sin(np.pi/2))
C_rr = (Srr - ne*Stt)/Ee
C_tt = (Stt - ne*Srr)/Ee
C_gr = 2*(1+ne)*Srt/Ee
print("C_rr   = %.12f" % C_rr)
print("C_tt   = %.12f" % C_tt)
print("C_grth = %.12f" % C_gr)
kappa = (3-ne)/(1+ne)
Cy = (KA/(2*G))*kappa/np.sqrt(2)
print("C_y90  = %.12f" % Cy)

ra = (KA/(10*sy))**2
M = max((1-ne)/Ee, 2*ne/Ee)
rb = (KA*M/0.05)**2
print("r_a = %.12f  r_b = %.12f  upper = min(ra,0.2)=%.12f" % (ra, rb, min(ra,0.2)))
print("sig_e at r: K_A/sqrt(r) ; ezz coeff = %.12f" % (-2*ne/Ee*KA))
