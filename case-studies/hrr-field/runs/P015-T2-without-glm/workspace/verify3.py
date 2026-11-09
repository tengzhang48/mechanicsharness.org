import numpy as np
from scipy.integrate import solve_bvp
import solve as S

n, alpha = 4.0, 2*np.pi/3
p0, q0 = 0.1561752603, -1.26307900
R0 = ((2-p0)+q0)/((2-p0)*(1-p0))

# vectorized f'''' RHS (same equation as shooting; solve_bvp does its own
# adaptive collocation + Newton on a global mesh -> independent discretization)
def vec_rhs(th, y, p):
    p = float(np.asarray(p).ravel()[0])
    c0 = p*(2.0-p); a = n*p*(2.0-n*p); m = n-1.0; c1 = 4.0*(1.0-n*p)
    f, f1, f2, f3 = y
    D = f2 + c0*f; Dp = f3 + c0*f1
    T = -(1.0-p)*f1; Tp = -(1.0-p)*f2; Tpp = -(1.0-p)*f3
    E2 = 0.75*D*D + 3.0*T*T
    E = np.sqrt(E2)
    Ep = (1.5*D*Dp + 6.0*T*Tp)/(2.0*E)
    num = m*((n-3.0)*Ep*Ep*D/E2 + D*(1.5*Dp*Dp + 6.0*Tp*Tp + 6.0*T*Tpp)/(2.0*E2)) \
          + 2.0*m*Ep*Dp/E + a*D - c1*(m*Ep*T/E + Tp)
    den = 1.0 + 0.75*m*D*D/E2
    Dpp = -num/den
    return np.vstack([f1, f2, f3, Dpp - c0*f2])

def bc(ya, yb, p):
    return np.array([ya[0]-1.0, ya[1], ya[3], yb[0], yb[1]])

# initial guess from shooting solution
sol0 = S.shoot(p0, q0, n, alpha, rtol=1e-12, atol=1e-14)
x_init = np.linspace(0, alpha, 400)
y_init = sol0.sol(x_init)

sol = solve_bvp(vec_rhs, bc, x_init, y_init, p=[p0], tol=1e-11, max_nodes=200000, verbose=1)
pb = float(sol.p[0])
qb = float(sol.sol(0.0)[2])
Rb = ((2-pb)+qb)/((2-pb)*(1-pb))
print(f"solve_bvp : p={pb:.10f} q={qb:.8f} R={Rb:.8f} status={sol.status} rms_max={sol.rms_residuals.max():.2e}")
print(f"vs shoot  : p={p0:.10f} q={q0:.8f} R={R0:.8f}")
print(f"diffs     : dp={abs(pb-p0):.2e}  dR={abs(Rb-R0):.2e}")

# also redo the crack check with solve_bvp as an extra independent confirmation
alpha_c = np.pi
p0c, q0c = 0.2, -0.65791
sol0c = S.shoot(p0c, q0c, 4.0, alpha_c, rtol=1e-12, atol=1e-14)
xc = np.linspace(0, alpha_c, 400)
yc = sol0c.sol(xc)
solc = solve_bvp(vec_rhs, lambda ya, yb, p: np.array([ya[0]-1.0, ya[1], ya[3], yb[0], yb[1]]),
                 xc, yc, p=[p0c], tol=1e-11, max_nodes=200000)
print(f"crack n=4 via solve_bvp: p={float(solc.p[0]):.10f} (target 0.2)")
