"""Independent check: treat the separable exponent t as an unknown nonlinear
eigenvalue and solve the two-point BVP  U'(0)=0, U(pi)=0, U(0)=1.
Must recover t = -1/(n+1) and the same angular shape as the IVP solution."""
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp

n = 3.5
a = (n - 1.0) / 2.0


def rhs_fields(t, U1, U2):
    c = t + 1.0
    k = (n * t + 1.0) * c
    Q = U1 * U1 + (U2 / c) ** 2
    U2p = -U1 * (k * Q + 2.0 * a * U2 * U2) / (Q + 2.0 * a * U2 * U2 / (c * c))
    return U2p


def fun(x, y, p):
    t = p[0]
    U2p = rhs_fields(t, y[0], y[1])
    return np.vstack((y[1], U2p))


def bc(ya, yb, p):
    return np.array([ya[1], yb[0], ya[0] - 1.0])


x = np.linspace(0.0, np.pi, 201)
yguess = np.vstack((np.cos(x / 2.0), -0.5 * np.sin(x / 2.0)))
sol = solve_bvp(fun, bc, x, yguess, p=[-0.2], tol=1e-10, max_nodes=200000)
t_bvp = sol.p[0]
print(f"BVP status: {sol.status} ({sol.message.strip()})")
print(f"BVP eigenvalue t = {t_bvp:.12f}   vs  -1/(n+1) = {-1.0/(n+1.0):.12f}"
      f"   diff = {abs(t_bvp + 1.0/(n+1.0)):.3e}")

# IVP reference (same equations, t fixed to -1/(n+1))
t_ref = -1.0 / (n + 1.0)


def ivp_rhs(th, y):
    return [y[1], rhs_fields(t_ref, y[0], y[1])]


ref = solve_ivp(ivp_rhs, (0.0, np.pi), [1.0, 0.0], dense_output=True,
                rtol=1e-12, atol=1e-14)
xs = np.linspace(0.0, np.pi, 501)
dU = np.max(np.abs(sol.sol(xs)[0] - ref.sol(xs)[0]))
dUp = np.max(np.abs(sol.sol(xs)[1] - ref.sol(xs)[1]))
print(f"shape comparison BVP vs IVP: max|dU| = {dU:.3e}, max|dU'| = {dUp:.3e}")
print(f"BVP residual max: {np.max(np.abs(sol.rms_residuals)) if hasattr(sol,'rms_residuals') else float('nan'):.3e}")
