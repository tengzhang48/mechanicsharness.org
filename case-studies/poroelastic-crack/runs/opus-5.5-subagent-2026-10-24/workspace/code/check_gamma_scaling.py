"""Check that the t=0+ tip factor depends only on gamma=(nu_u-nu)/(1-nu):
W(s) = s*K1 + 2*beta*s*K2 (loading sigma_n=-1, p=p_u on faces) at large s for several
parameter sets with gamma = 0.25."""
import numpy as np
from crack_solver2 import Params, solve
sets = [(0.2, 0.4, 0.89), (0.2, 0.4, 0.5), (0.0, 0.25, 0.7), (0.1, 0.325, 1.0)]
for (nu, nuu, al) in sets:
    P = Params(nu, nuu, al)
    out = []
    for s in [1e2, 1e4]:
        o = solve(s, 40, P)
        W = s*o['L1']['K'] + 2*P.beta*s*o['L2']['K']
        VW = s*o['L1']['V'] + 2*P.beta*s*o['L2']['V']
        out.append((s, W, (VW - np.pi*(1 - nuu))*np.sqrt(s)/(1 - nuu)))
    print(f"nu={nu}, nu_u={nuu}, alpha={al}, gamma={P.gam:.4f}: " + "  ".join(f"s={s:.0e}: W={W:.7f}, (VW-Vu)sqrt(s)/(1-nu_u)={c:.6f}" for s, W, c in out))
