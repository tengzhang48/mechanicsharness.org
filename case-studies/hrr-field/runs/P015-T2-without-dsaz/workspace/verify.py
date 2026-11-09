import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

def check(n, alpha, s0, M=500):
    sol = d.run(n, alpha, s0, M)
    s = sol.p[0]
    th = sol.x
    U, Up, V, P = sol.y
    # verify the constitutive/equilibrium eqs pointwise using our own eqs()
    # recompute V' from inc and compare to stored V' via finite diff
    Vp_inc = -(2 - s*n)*U
    Vp_num = np.gradient(V, th)
    # residual of equilibrium eqs using solve_Upp_Pp on the stored field
    Upp, Pp = d.solve_Upp_Pp(U, Up, V, P, s, n)
    Upp_num = np.gradient(Up, th)
    Pp_num = np.gradient(P, th)
    print(f"n={n} alpha={alpha:.4f} s0={s0} -> s={s:.6f} status={sol.status} rms={np.max(sol.rms_residuals):.3e}")
    print(f"   inc err max|Vp_inc - Vp_num| = {np.max(np.abs(Vp_inc - Vp_num)):.3e}")
    print(f"   U'' err max = {np.max(np.abs(Upp - Upp_num)):.3e}")
    print(f"   P'  err max = {np.max(np.abs(Pp - Pp_num)):.3e}")
    Srr, Stt, Srt = d.stresses(U, Up, V, P, s, n)
    print(f"   Srr(0)={Srr[0]:+.5f} Stt(0)={Stt[0]:+.5f} R={Srr[0]/Stt[0]:+.5f}")
    print(f"   Stt(alpha)={Stt[-1]:.2e} Srt(alpha)={Srt[-1]:.2e}")

def check_sol(sol, n, alpha, label):
    s = sol.p[0]
    th = sol.x
    U, Up, V, P = sol.y
    Vp_inc = -(2 - s*n)*U
    Vp_num = np.gradient(V, th)
    Upp, Pp = d.solve_Upp_Pp(U, Up, V, P, s, n)
    Upp_num = np.gradient(Up, th)
    Pp_num = np.gradient(P, th)
    print(f"[{label}] n={n} alpha={alpha:.4f} -> s={s:.6f} status={sol.status} rms={np.max(sol.rms_residuals):.3e}")
    print(f"   max|Vp_inc-Vp_num|={np.max(np.abs(Vp_inc-Vp_num)):.3e}  max|U''err|={np.max(np.abs(Upp-Upp_num)):.3e}  max|P'err|={np.max(np.abs(Pp-Pp_num)):.3e}")
    Srr, Stt, Srt = d.stresses(U, Up, V, P, s, n)
    print(f"   Srr(0)={Srr[0]:+.5f} Stt(0)={Stt[0]:+.5f} R={Srr[0]/Stt[0]:+.5f}  Stt(a)={Stt[-1]:.2e} Srt(a)={Srt[-1]:.2e}")

if __name__ == "__main__":
    alpha = 2*np.pi/3; M = 500
    # build physical (Mode A) via continuation from n=1
    prev = d.run(1.0, alpha, 0.3843, M)
    for n in np.arange(1.25, 4.0001, 0.25):
        thn = np.linspace(0, alpha, M)
        Y = np.zeros((4, M))
        for k in range(4): Y[k] = np.interp(thn, prev.x, prev.y[k])
        bcw = lambda ya, yb, p, n=n: d.bc(ya, yb, p[0], n)
        rhsw = lambda x, y, p, n=n: d.rhs(x, y, p[0], n)
        prev = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([prev.p[0]]), tol=1e-10, max_nodes=10000)
    check_sol(prev, 4.0, alpha, "MODE A physical")
    solB = d.run(4.0, alpha, 0.35, M)
    check_sol(solB, 4.0, alpha, "MODE B other")