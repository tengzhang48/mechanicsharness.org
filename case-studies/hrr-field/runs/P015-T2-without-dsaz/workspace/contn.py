import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

alpha = 2*np.pi/3
n_list = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
prev = None
s_guess = 0.3843
for n in n_list:
    # seed from previous solution if available
    if prev is None:
        sol = d.run(n, alpha, s_guess, 400)
    else:
        # build guess on same grid by interpolation
        from numpy import interp
        th_new = np.linspace(0, alpha, 400)
        Y = np.zeros((4, 400))
        for k in range(4):
            Y[k] = interp(th_new, prev.x, prev.y[k])
        a0 = np.array([prev.p[0]])
        def bcw(ya, yb, p, n=n):
            return d.bc(ya, yb, p[0], n)
        def rhsw(x, y, p, n=n):
            return d.rhs(x, y, p[0], n)
        sol = d.solve_bvp(rhsw, bcw, th_new, Y, p=a0, tol=1e-9, max_nodes=6000)
    s = sol.p[0]
    Srr, Stt, Srt = d.stresses(sol.y[0], sol.y[1], sol.y[2], sol.y[3], s, n)
    R = Srr[0]/Stt[0]
    print(f"n={n:4.2f}  s={s:.6f}  status={sol.status}  rms={np.max(sol.rms_residuals):.2e}  R=Srr(0)/Stt(0)={R:+.5f}")
    prev = sol