"""Final assembly: Part A and Part B from the verified first-order operator.

First-order operator (global-basis, per z-Fourier mode n, sigma = sign(n)):
  dK1^g        = |k| da_hat ( -1/2 ) K1^0
  dK2^g        = |k| da_hat [ (-1/2+g)K2^0 + i sigma g K3^0 ]
  dK3^g        = |k| da_hat [ (-1/2-g)K3^0 + i sigma g K2^0 ]
  g = nu/(2-nu)
Problem-basis K's (local n,t):  K2^prob = K2^g - a' K3^0 ; K3^prob = K3^g + a' K2^0
with a' = i k da (mode-wise i sigma |k| da_hat).
Basis rotation validated by the exact inclined-front anti-plane test.

Part A: da = A cos(kz), k = 2 pi / lam.
Part B: per-mode equilibrium G = Gamma:
  kappa1 = -pi|n| f_n
  kappa2 = 2 pi|n| T2(sigma) f_n,  T2(+) = (d2)(0.6) + i g(0.8) - i(0.8)
  kappa3 = 2 pi|n| T3(sigma) f_n,  T3(+) = (d3)(0.8) + i g(0.6) + i(0.6)
  2[w1(k1 + 0.6 k2) + w3 0.8 k3]/S0 = ghat_n  ->  f_n.
"""
import numpy as np
import json

NU = 0.3
W1 = 1 - NU**2
W3 = 1 + NU
R2 = 0.6; R3 = 0.8
S0 = W1*(1 + R2**2) + W3*R3**2
GAM = NU/(2 - NU)
D2 = -0.5 + GAM
D3 = -0.5 - GAM

T2p = D2*R2 + 1j*GAM*R3 - 1j*R3      # kappa2 per (2 pi n f_n), n>0
T3p = D3*R3 + 1j*GAM*R2 + 1j*R2      # kappa3 per (2 pi n f_n), n>0

# ---------------- Part A ----------------
C1 = -np.pi
S1 = 0.0
C2 = 2*np.pi*T2p.real
S2 = -2*np.pi*T2p.imag
C3 = 2*np.pi*T3p.real
S3 = -2*np.pi*T3p.imag
CG = 2*(W1*(C1 + R2*C2) + W3*R3*C3)/S0
SG = 2*(W1*R2*S2 + W3*R3*S3)/S0
print("Part A:")
for nm, v in [('C1', C1), ('S1', S1), ('C2', C2), ('S2', S2), ('C3', C3), ('S3', S3),
              ('CG', CG), ('SG', SG)]:
    print(f"  {nm} = {v:.8f}")

# ---------------- Part B ----------------
K10 = np.sqrt(W1/S0)
print(f"\nPart B: K10_scaled = {K10:.10f}")

def Bcoef(n):
    """operator bracket per mode n: kappa vector = 2 pi |n| f_n * (r_y, T2, T3)"""
    sg = np.sign(n)
    T2 = D2*R2 + 1j*sg*GAM*R3 - 1j*sg*R3
    T3 = D3*R3 + 1j*sg*GAM*R2 + 1j*sg*R2
    return (-0.5, T2, T3)

q = 0.5
NMAX = 60
nvals = np.arange(1, NMAX+1)
fhat = {}
for n in nvals:
    gh = q**n
    ry, T2, T3 = Bcoef(n)
    B = W1*(ry + R2*T2) + W3*R3*T3
    fhat[n] = gh*S0/(4*np.pi*n*B)
    fhat[-n] = np.conj(fhat[n])

def series_vec(coef_by_n, svals):
    svals = np.asarray(svals)
    tot = np.zeros(svals.shape, complex)
    for n in nvals:
        ph = np.exp(2j*np.pi*n*svals)
        tot += coef_by_n(n)*ph + np.conj(coef_by_n(n))*np.conj(ph)
    return tot

# check zero mean of f
print("  f mean check (should be ~0):", series_vec(lambda n: fhat[n], np.linspace(0, 1, 200, endpoint=False)).mean())

f_of_s  = lambda s: series_vec(lambda n: fhat[n], s).real
k1_of_s = lambda s: series_vec(lambda n: -np.pi*n*fhat[n], s).real
k2_of_s = lambda s: series_vec(lambda n: 2*np.pi*n*Bcoef(n)[1]*fhat[n], s).real
k3_of_s = lambda s: series_vec(lambda n: 2*np.pi*n*Bcoef(n)[2]*fhat[n], s).real

sk = np.arange(16)/16.0
f16 = [f_of_s(s) for s in sk]
k16_1 = [k1_of_s(s) for s in sk]
k16_2 = [k2_of_s(s) for s in sk]
k16_3 = [k3_of_s(s) for s in sk]
print("  f(k/16):", [f"{v:+.6f}" for v in f16])
print("  k1:", [f"{v:+.6f}" for v in k16_1])
print("  k2:", [f"{v:+.6f}" for v in k16_2])
print("  k3:", [f"{v:+.6f}" for v in k16_3])

sf = np.linspace(-0.5, 0.5, 200001)
k3v = k3_of_s(sf)
fv = f_of_s(sf)
i3 = int(np.argmax(k3v))
print(f"  kappa3_max = {k3v[i3]:.8f} at s = {sf[i3]:.6f}")
print(f"  f max = {fv.max():.8f} at s={sf[np.argmax(fv)]:.4f}, f min = {fv.min():.8f}")
print(f"  front_peak_to_peak = {fv.max()-fv.min():.8f}")

# pointwise equilibrium check: 2[w1(k1+0.6k2)+w3 0.8k3]/S0 must equal g(s)
gs = lambda s: (1-q**2)/(1 - 2*q*np.cos(2*np.pi*s) + q**2) - 1.0
ss = np.array([0.0, 0.03, 0.11, 0.25, 0.4, 0.5])
lhs = 2*(W1*(k1_of_s(ss) + R2*k2_of_s(ss)) + W3*R3*k3_of_s(ss))/S0
rhs = gs(ss)
print("  equilibrium check lhs-g(s):", np.array_str(lhs-rhs, precision=8))

results = {
  "partA": {"C1": float(C1), "S1": float(S1), "C2": float(C2), "S2": float(S2),
            "C3": float(C3), "S3": float(S3), "CG": float(CG), "SG": float(SG)},
  "partB": {"K10_scaled": float(K10), "front_peak_to_peak": float(fv.max()-fv.min()),
            "f": [float(v) for v in f16], "kappa1": [float(v) for v in k16_1],
            "kappa2": [float(v) for v in k16_2], "kappa3": [float(v) for v in k16_3],
            "kappa3_max": float(k3v[i3]), "s_kappa3_max": float(sf[i3])}
}
with open('/work/results.json', 'w') as fh:
    json.dump(results, fh, indent=1)
print("\nsaved /work/results.json")
