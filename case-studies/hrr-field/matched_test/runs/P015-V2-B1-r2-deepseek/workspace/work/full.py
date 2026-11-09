import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

E = 70000.0
nu = 0.30
sy = 10.0
J = 0.800
Czz0 = -0.005000

invE = 1.0 / E
A_hyd = (1 - 2 * nu) / (3 * E)

def A_dev(T):
    return 1.5 / T + (nu - 0.5) / E

def eff_props(T):
    G = 1.0 / (2.0 * A_dev(T))
    K = 1.0 / (9.0 * A_hyd)
    Eeff = 9.0 * K * G / (3.0 * K + G)
    nueff = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return G, K, Eeff, nueff

def Czz_of_T(T):
    G, K, Eeff, nueff = eff_props(T)
    KA = np.sqrt(J * Eeff / (2.0 * np.pi))
    return -2.0 * nueff / Eeff * KA

T = brentq(lambda T: Czz_of_T(T) - Czz0, 1000.0, 60000.0, xtol=1e-13, rtol=1e-15)
G, K, Eeff, nueff = eff_props(T)
KA = np.sqrt(J * Eeff / (2.0 * np.pi))
print("T          = %.9f MPa" % T)
print("E_eff      = %.9f MPa  (equals T ? %s)" % (Eeff, abs(Eeff - T) < 1e-6))
print("nu_eff     = %.9f" % nueff)
print("G_eff      = %.9f" % G)
print("K_A        = %.9f MPa*mm^0.5" % KA)
print("nu_eff pred= %.9f" % (0.5 + T * (nu - 0.5) / E))

# ---- angular stress functions (sigma_ij = r^-1/2 * Sigma_ij(theta)), mode I, plane stress
th = np.linspace(-np.pi, np.pi, 400001)
c2, c32 = np.cos(th * 0.5), np.cos(1.5 * th)
s2, s32 = np.sin(th * 0.5), np.sin(1.5 * th)
Srr = KA * (1.25 * c2 - 0.25 * c32)
Stt = KA * (0.75 * c2 + 0.25 * c32)
Srt = KA * (0.25 * s2 + 0.25 * s32)

def strains_angular(theta):
    c2, c32 = np.cos(theta * 0.5), np.cos(1.5 * theta)
    s2, s32 = np.sin(theta * 0.5), np.sin(1.5 * theta)
    srr = (1.25 * c2 - 0.25 * c32)
    stt = (0.75 * c2 + 0.25 * c32)
    srt = (0.25 * s2 + 0.25 * s32)
    S = srr + stt
    err = A_hyd * S + A_dev(T) * (srr - S / 3.0)
    ett = A_hyd * S + A_dev(T) * (stt - S / 3.0)
    ert = A_dev(T) * srt
    ezz = A_hyd * S + A_dev(T) * (-S / 3.0)
    return srr, stt, srt, err, ett, ert, ezz

# ---- part 2 : theta = 60 deg
t60 = np.pi / 3.0
srr, stt, srt, err, ett, ert, ezz = strains_angular(np.array([t60]))
print("\nPart 2 at theta=60 deg:")
print("  Sigma_rr/K_A=%.9f Sigma_tt/K_A=%.9f Sigma_rt/K_A=%.9f" % (srr[0], stt[0], srt[0]))
C_rr = KA * err[0]
C_tt = KA * ett[0]
C_grg = KA * 2.0 * ert[0]
print("  C_rr  = %.9f mm^1/2" % C_rr)
print("  C_tt  = %.9f mm^1/2" % C_tt)
print("  C_grth= %.9f mm^1/2" % C_grg)
print("  (alt via E_eff,nu_eff: C_rr=%.9f C_tt=%.9f C_grth=%.9f)"
      % (KA / Eeff * (srr[0] - nueff * stt[0]),
         KA / Eeff * (stt[0] - nueff * srr[0]),
         KA * (1 + nueff) / Eeff * srt[0] * 2))
print("  C_zz(0) check = %.9f" % (KA * ezz[0] / 1.0 * 0 + KA * (A_hyd - A_dev(T) / 3.0) * 2.0))
