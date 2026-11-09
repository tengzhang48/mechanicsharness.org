import numpy as np
from scipy.integrate import solve_bvp
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

def report(n, alpha, s0):
    sol = d.run(n, alpha, s0, 300)
    s = sol.p[0]
    th = sol.x
    U = sol.y[0]; Up = sol.y[1]; V = sol.y[2]; P = sol.y[3]
    Srr, Stt, Srt = d.stresses(U, Up, V, P, s, n)
    print(f"--- n={n} alpha={alpha:.4f} s0={s0}: s={s:.5f} status={sol.status} nmax={sol.x.size}")
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        i = int(frac*(len(th)-1))
        print(f"   th/a={frac:.2f} U={U[i]:+.4f} V={V[i]:+.4f} Srr={Srr[i]:+.4f} Stt={Stt[i]:+.4f} Srt={Srt[i]:+.4f} P={P[i]:+.4f}")
    # sigma_e (proportional to e2^{1/2}) and ratio at 0
    print(f"   Srr(0)/Stt(0)={Srr[0]/(Stt[0]+1e-30):+.4f}")
    return s, (U, V, P, Srr, Stt, Srt, th)

if __name__ == "__main__":
    report(4.0, 2*np.pi/3, 0.20)
    report(4.0, 2*np.pi/3, 0.35)