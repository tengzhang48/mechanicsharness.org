import numpy as np
import json

nu = 0.3
q = 0.5
K1, K2, K3 = 1.0, 0.6, 0.8
w = 1 - nu**2
v = 1 + nu
G0 = w*(K1**2 + K2**2) + v*K3**2

# Rice (1985) first-order front-perturbation response:
#   delta K_a = K_b^0 [ -A_ab |k| a + B_ab a' ]
# mode I stable, modes II & III unstable; tilt couples II<->III
a1, a2, a3 = 0.5, -0.5, -0.5
B23, B32 = -1.0, 1.0

Cadv = 2*(w*(a1*K1**2 + a2*K2**2) + v*a3*K3**2)/G0
Ctilt = 2*K2*K3*(w*B23 + v*B32)/G0
Delta = Cadv**2 + Ctilt**2

# g(s) = (1-q^2)/(1-2q cos 2pi s + q^2) - 1 = sum_{n>=1} 2 q^n cos(2 pi n s)
N = 4000
n = np.arange(1, N+1, dtype=float)
qp = q**n
fc = (-Cadv)*qp/(np.pi*n*Delta)
fs = ( Ctilt)*qp/(np.pi*n*Delta)

def fval(s):
    th = 2*np.pi*s
    return np.sum(fc*np.cos(n*th) + fs*np.sin(n*th))
def fpval(s):
    th = 2*np.pi*s
    return np.sum(-2*np.pi*n*fc*np.sin(n*th) + 2*np.pi*n*fs*np.cos(n*th))
def Hval(s):
    th = 2*np.pi*s
    return np.sum(2*np.pi*n*fc*np.cos(n*th) + 2*np.pi*n*fs*np.sin(n*th))
def gval(s):
    th = 2*np.pi*s
    return np.sum(2*qp*np.cos(n*th))

def kappas(s):
    H = Hval(s); Fp = fpval(s)
    k1 = -a1*(K1/K1)*H
    k2 = -a2*(K2/K1)*H + B23*(K3/K1)*Fp
    k3 = -a3*(K3/K1)*H + B32*(K2/K1)*Fp
    return k1,k2,k3

s16 = np.arange(16)/16.0
f16 = np.array([fval(s) for s in s16])
k1_16 = np.empty(16); k2_16 = np.empty(16); k3_16 = np.empty(16)
for i,s in enumerate(s16):
    k1_16[i], k2_16[i], k3_16[i] = kappas(s)

# fine grid for front peak-to-peak and kappa3 max
sf = np.linspace(-0.5, 0.5, 200000, endpoint=False)
ff = np.array([fval(s) for s in sf])
k3f = np.array([kappas(s)[2] for s in sf])
imax = int(np.argmax(k3f))
front_pp = float(ff.max() - ff.min())
k3max = float(k3f[imax])
smax = float(sf[imax])

# part A coefficients
C1 = -2*np.pi*a1*(K1/K1)
S1 = 0.0
C2 = -2*np.pi*a2*(K2/K1)
S2 = -2*np.pi*B23*(K3/K1)
C3 = -2*np.pi*a3*(K3/K1)
S3 = -2*np.pi*B32*(K2/K1)
CG = -2*np.pi*Cadv
SG = -2*np.pi*Ctilt
K10_scaled = float(np.sqrt((1-nu**2)/G0))

# consistency check: G from kappas equals g
print("G consistency (should be ~0):")
for i,s in enumerate(s16):
    k1,k2,k3 = kappas(s)
    dG = (2*w*K1*k1 + 2*w*K2*k2 + 2*v*K3*k3)/G0
    print(f"  s={s:.4f} dG={dG:.10f} g={gval(s):.10f} diff={dG-gval(s):.2e}")
print("Cadv,Ctilt,Delta",Cadv,Ctilt,Delta)
print("f16",repr(f16))
print("k1",repr(k1_16))
print("k2",repr(k2_16))
print("k3",repr(k3_16))
print("front_pp",front_pp,"k3max",k3max,"smax",smax)
print("PartA",C1,S1,C2,S2,C3,S3,CG,SG)
print("K10",K10_scaled)

results = {
  "partA": {"C1": float(C1), "S1": float(S1), "C2": float(C2), "S2": float(S2),
            "C3": float(C3), "S3": float(S3), "CG": float(CG), "SG": float(SG)},
  "partB": {"K10_scaled": K10_scaled, "front_peak_to_peak": front_pp,
            "f": [float(x) for x in f16],
            "kappa1": [float(x) for x in k1_16],
            "kappa2": [float(x) for x in k2_16],
            "kappa3": [float(x) for x in k3_16],
            "kappa3_max": k3max, "s_kappa3_max": smax}
}
with open("results.json","w") as fh:
    json.dump(results, fh, indent=2)
print("wrote results.json")
