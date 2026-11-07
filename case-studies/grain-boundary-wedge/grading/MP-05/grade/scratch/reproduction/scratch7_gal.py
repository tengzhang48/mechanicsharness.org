import numpy as np
from solver3 import galerkin_eigs
from core import kernel_matrix, f_matrix

for N in [16, 24, 32, 40, 48]:
    mu, Z, M, G, X = galerkin_eigs(N)
    sel = [(m.real, i) for i, m in enumerate(mu)
           if abs(m.imag) < 1e-8*max(1.0, abs(m.real)) and m.real < 0]
    sel.sort(key=lambda t: -t[0])   # largest (least negative) = slowest decay
    line = f"N={N:3d}: {len(sel)} real decaying; slowest mu: "
    line += "  ".join(f"{m:.6f}" for m, i in sel[:4])
    print(line)
