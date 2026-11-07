import numpy as np
from solver4 import build, spectrum
from core import f_matrix

for N in [16, 24, 32, 40]:
    E, c0, B, xq = build(N)
    mu, Z = spectrum(E)
    sel = sorted([m.real for m in mu if abs(m.imag) < 1e-7*max(1,abs(m.real)) and m.real < -1e-6])
    sel.reverse()
    print(f"N={N:3d}: slowest real mu: {['%.6f'%v for v in sel[:5]]}")

# trajectory at x=0.5
from scipy.linalg import expm
E, c0, B, xq = build(32)
bo = f_matrix(np.array([0.5]), 32)[0]
ss = np.linspace(0.0, 0.02, 41)
vals = [bo @ (expm(E*s) @ c0) for s in ss]
print("sigma(0.5,s) sample:", [f"{v:.6f}" for v in vals[::8]])
