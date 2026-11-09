import numpy as np
from scipy.integrate import solve_bvp

# Displacement formulation of the self-similar incompressible power-law field.
#   u_r = r^{1-sn} U(theta), u_theta = r^{1-sn} V(theta),  pressure r^{-sn} P(theta)
# Incompressibility (shapes): (2-sn)U + V' = 0
# Deviatoric strain shapes:
#   d_rr = (1-sn)U
#   d_tt = U + V'          (= -(1-sn)U by inc)
#   d_rt = (1/2)(U' - sn V)
# Constitutive (B=1):  K = e2^{(1-n)/(2n)},  s_ij = (2/3) K d_ij,
#   e2 = (2/3)( d_rr^2 + d_tt^2 + 2 d_rt^2 )
# Stresses (shape):  S_rr = s_rr + P, S_tt = s_tt + P, S_rt = s_rt
# Equilibrium (r-,theta-):
#   eqr: -sn S_rr + (S_rr - S_tt) + S_rt' = 0
#   eqt: (-sn+2) S_rt + S_tt' = 0
# We use inc to replace V' = -(2-sn)U and V'' = -(2-sn)U'.
# Then (eqr,eqt) are LINEAR in the unknowns (U'', P'); solve the 2x2.

def stresses(U, Up, V, P, s, n):
    d_rr = (1 - s*n)*U
    d_tt = -(1 - s*n)*U        # by incompressibility
    d_rt = 0.5*(Up - s*n*V)
    e2 = (2.0/3.0)*(d_rr**2 + d_tt**2 + 2*d_rt**2) + 1e-300
    K = e2**((1-n)/(2.0*n))
    s_rr = (2.0/3.0)*K*d_rr
    s_tt = (2.0/3.0)*K*d_tt
    s_rt = (2.0/3.0)*K*d_rt
    return s_rr + P, s_tt + P, s_rt

def eqs(U, Up, Upp, V, P, Pp, s, n):
    # evaluate (eqr, eqt) given U'' and P'
    S_rr, S_tt, S_rt = stresses(U, Up, V, P, s, n)
    # S_rt' and S_tt' : need derivatives. Build with small central differences
    # by perturbing theta implicitly. Use analytic: d/dtheta of the shapes.
    h = 1e-6
    # helper to evaluate stress shapes at (U+..., V+..., P+...)
    def sh(Ux, Upx, Vx):
        d_rr = (1 - s*n)*Ux
        d_tt = -(1 - s*n)*Ux
        d_rt = 0.5*(Upx - s*n*Vx)
        e2 = (2.0/3.0)*(d_rr**2 + d_tt**2 + 2*d_rt**2) + 1e-300
        K = e2**((1-n)/(2.0*n))
        return (2.0/3.0)*K*d_rr, (2.0/3.0)*K*d_tt, (2.0/3.0)*K*d_rt
    # forward difference in theta: U->U+h*Up+h^2/2 Upp ; Up->Up+h*Upp ; V->V+h*V' ; V'=-(2-sn)U
    def dtheta():
        U1 = U + h*Up + 0.5*h*h*Upp
        Up1 = Up + h*Upp
        V1 = V + h*(-(2 - s*n)*U)
        P1 = P + h*Pp
        s_rr1, s_tt1, s_rt1 = sh(U1, Up1, V1)
        Srt1 = s_rt1; Stt1 = s_tt1 + P1
        return Srt1, Stt1
    S_rt1, S_tt1 = dtheta()
    S_rt_p = (S_rt1 - S_rt)/h
    S_tt_p = (S_tt1 - S_tt)/h
    eqr = -s*n*S_rr + (S_rr - S_tt) + S_rt_p
    eqt = (-s*n + 2)*S_rt + S_tt_p
    return eqr, eqt

def solve_Upp_Pp(U, Up, V, P, s, n):
    # solve linear 2x2 for (Upp, Pp) using 3 evaluations (linear system),
    # vectorized over the node axis.
    U, Up, V, P = np.atleast_1d(U), np.atleast_1d(Up), np.atleast_1d(V), np.atleast_1d(P)
    z = np.zeros_like(U)
    f0 = eqs(U, Up, z, V, P, z, s, n)
    fU = eqs(U, Up, z + 1.0, V, P, z, s, n)
    fP = eqs(U, Up, z, V, P, z + 1.0, s, n)
    a11 = fU[0] - f0[0]; a12 = fP[0] - f0[0]
    a21 = fU[1] - f0[1]; a22 = fP[1] - f0[1]
    det = a11*a22 - a12*a21
    det = np.where(np.abs(det) < 1e-30, 1e-30, det)
    b1 = -f0[0]; b2 = -f0[1]
    Upp = (b1*a22 - b2*a12)/det
    Pp = (a11*b2 - a21*b1)/det
    return Upp, Pp

def rhs(th, Y, s, n):
    # Y = [U, Up, V, P]
    U, Up, V, P = Y
    Upp, Pp = solve_Upp_Pp(U, Up, V, P, s, n)
    Vp = -(2 - s*n)*U
    return np.vstack([Up, Upp, Vp, Pp])

def bc(ya, yb, s, n):
    Sa = stresses(ya[0], ya[1], ya[2], ya[3], s, n)
    Sb = stresses(yb[0], yb[1], yb[2], yb[3], s, n)
    return np.array([
        ya[1],        # U'(0)=0   (mode-I evenness)
        ya[2],        # V(0)=0    (mode-I oddness)
        Sb[2],        # S_rt(alpha)=0  traction free
        Sb[1],        # S_tt(alpha)=0  traction free
        ya[0] - 1.0,  # amplitude gauge
    ])

def run(n, alpha, s0, M=200):
    th = np.linspace(0, alpha, M)
    Y = np.zeros((4, M))
    Y[0] = np.cos(np.pi*th/(2*alpha))
    Y[2] = -np.sin(np.pi*th/(2*alpha))
    a0 = np.array([s0])
    def bcw(ya, yb, p): return bc(ya, yb, p[0], n)
    def rhsw(x, y, p): return rhs(x, y, p[0], n)
    sol = solve_bvp(rhsw, bcw, th, Y, p=a0, tol=1e-8, max_nodes=4000)
    return sol

if __name__ == "__main__":
    import sys
    n = float(sys.argv[1]); alpha = float(sys.argv[2]); s0 = float(sys.argv[3])
    M = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    sol = run(n, alpha, s0, M)
    print("n=%g alpha=%.4f s0=%.3f -> status=%s s=%.6f rmsres=%.2e"
          % (n, alpha, s0, sol.status, sol.p[0], np.max(sol.rms_residuals)))
    print("U(0)=%.4f U(a)=%.4f nodes=%d" % (sol.y[0][0], sol.y[0][-1], sol.x.size))