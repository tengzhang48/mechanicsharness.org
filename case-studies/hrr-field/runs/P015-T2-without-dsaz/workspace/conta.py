import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

def cont(n, a_list, s0, M=400):
    prev = None
    for alpha in a_list:
        if prev is None:
            sol = d.run(n, alpha, s0, M)
        else:
            th_new = np.linspace(0, alpha, M)
            Y = np.zeros((4, M))
            for k in range(4):
                Y[k] = np.interp(th_new, prev.x, prev.y[k])
            a0 = np.array([prev.p[0]])
            bcw = lambda ya, yb, p: d.bc(ya, yb, p[0], n)
            rhsw = lambda x, y, p: d.rhs(x, y, p[0], n)
            sol = d.solve_bvp(rhsw, bcw, th_new, Y, p=a0, tol=1e-9, max_nodes=8000)
        s = sol.p[0]
        Srr, Stt, Srt = d.stresses(sol.y[0], sol.y[1], sol.y[2], sol.y[3], s, n)
        R = Srr[0]/Stt[0]
        print(f"n={n} alpha={alpha:.5f} s={s:.6f} R={R:+.5f} status={sol.status} rms={np.max(sol.rms_residuals):.2e}")
        prev = sol
    return prev

if __name__ == "__main__":
    import sys
    n = float(sys.argv[1])
    a0 = 2*np.pi/3; a1 = np.pi*0.998
    a_list = list(np.linspace(a0, a1, 14))
    cont(n, a_list, 0.0876 if n > 3 else 0.3843)