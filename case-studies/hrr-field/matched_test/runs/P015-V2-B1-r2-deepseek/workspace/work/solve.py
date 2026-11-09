import numpy as np
from scipy.optimize import brentq

mm = 1.0
N = 1.0
MPa = 1.0

E = 70000.0
nu = 0.30
sy = 10.0
J = 0.800
Czz0 = -0.005000

invE = 1.0 / E
# compliance split:  eps_ij = A_hyd * skk * dij + A_dev * s_ij
A_hyd = (1 - 2 * nu) / (3 * E)
# A_dev = 1.5/T + (nu - 0.5)/E
def A_dev(T):
    return 1.5 / T + (nu - 0.5) / E

def eff_props(T):
    G = 1.0 / (2.0 * A_dev(T))          # shear modulus
    K = 1.0 / (9.0 * A_hyd)             # bulk modulus (since eps_kk = 3 A_hyd skk, skk=3 sig_h)
    Eeff = 9.0 * K * G / (3.0 * K + G)
    nueff = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return G, K, Eeff, nueff

# C_zz(0) = -2 * nu_eff/E_eff * K_A  ;  K_A = sqrt(J*E_eff/(2*pi))
def Czz_of_T(T):
    G, K, Eeff, nueff = eff_props(T)
    KA = np.sqrt(J * Eeff / (2.0 * np.pi))
    return -2.0 * nueff / Eeff * KA

def f(T):
    return Czz_of_T(T) - Czz0

print("Czz for a few T:")
for T in [10000, 30000, 50000, 60000, 69000]:
    print(T, Czz_of_T(T), eff_props(T))

Ts = np.linspace(1.0, E - 1e-6, 200000)
vals = np.array([f(T) for T in Ts])
# find bracket
br = None
for i in range(len(Ts) - 1):
    if vals[i] * vals[i+1] < 0:
        br = (Ts[i], Ts[i+1])
        break
print("bracket", br)
T = brentq(f, *br, xtol=1e-12, rtol=1e-14)
print("T = %.10f MPa" % T)
G, K, Eeff, nueff = eff_props(T)
KA = np.sqrt(J * Eeff / (2 * np.pi))
print("G_eff=%.6f K=%.6f E_eff=%.6f nu_eff=%.8f K_A=%.10f" % (G, K, Eeff, nueff, KA))
print("check Czz =", -2*nueff/Eeff*KA, "J ->", 2*np.pi*KA**2/Eeff)
print("A_dev=%.6e A_hyd=%.6e  P=A_hyd-A_dev/3=%.6e" % (A_dev(T), A_hyd, A_hyd - A_dev(T)/3))
print("2KA*(Ahyd-Adev/3) =", 2*KA*(A_hyd - A_dev(T)/3))
