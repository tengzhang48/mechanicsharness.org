import numpy as np
from scipy.integrate import quad
from scipy.integrate import solve_bvp
import sympy as sp

nu = 0.3
Ebar = 1.0
s0 = 1.0
ell = 1.0
h = 1.0
k = h*Ebar/ell**2      # k/(h Ebar) = 1/ell^2

def u2_exact(x, H):
    B = -s0*ell/(Ebar*np.cosh(H/(2*ell)))
    return B*np.sinh((x - H/2)/ell)

def u2p_exact(x, H):
    B = -s0*ell/(Ebar*np.cosh(H/(2*ell)))
    return B/ell*np.cosh((x - H/2)/ell)

def W_closed(H):
    return -h*s0**2*ell/Ebar*np.tanh(H/(2*ell))

def W_quad(H):
    f = lambda x: h*s0*u2p_exact(x,H) + 0.5*h*Ebar*u2p_exact(x,H)**2 + 0.5*k*u2_exact(x,H)**2
    val,_ = quad(f, 0, H, limit=200)
    return val

print("=== CHECK 1: energy integral by direct quadrature vs closed form ===")
for H in [0.1,0.25,0.5,1.0,2.0,5.0]:
    wq, wc = W_quad(H), W_closed(H)
    rel = abs(wq-wc)/abs(wc)
    print(f"H={H:5.2f}  W_quad={wq: .10e}  W_closed={wc: .10e}  rel.err={rel:.2e}  tol1e-9 -> {'PASS' if rel<1e-9 else 'FAIL'}")

print()
print("=== CHECK 2: wake ODE residual u2'' = u2/ell^2 (finite difference) ===")
H=1.0
for xtest in [0.1,0.5,0.9]:
    eps=1e-5
    u2pp=(u2_exact(xtest+eps,H)-2*u2_exact(xtest,H)+u2_exact(xtest-eps,H))/eps**2
    res=abs(u2pp-u2_exact(xtest,H)/ell**2)
    print(f"x={xtest}  resid={res:.2e}  tol1e-4 -> {'PASS' if res<1e-4 else 'FAIL'}")

print()
print("=== CHECK 3: traction-free strip-edge conditions ===")
for H in [0.3,1.0,3.0]:
    s22_0 = s0 + Ebar*u2p_exact(0.0,H)
    s22_H = s0 + Ebar*u2p_exact(H,H)
    print(f"H={H}: sigma22(0)={s22_0:.2e}, sigma22(H)={s22_H:.2e}, tol1e-10 -> {'PASS' if max(abs(s22_0),abs(s22_H))<1e-10 else 'FAIL'}")

print()
print("=== CHECK 4: numerical BVP solve (solve_bvp) vs closed form ===")
def ode(x,y,H):
    return np.vstack((y[1], y[0]/ell**2))
def bc(ya,yb,H):
    return np.array([ya[1]+s0/Ebar, yb[1]+s0/Ebar])
for H in [0.5,1.0,2.0]:
    xg=np.linspace(0,H,101)
    y0=np.zeros((2,xg.size))
    sol=solve_bvp(lambda x,y: ode(x,y,H), lambda ya,yb: bc(ya,yb,H), xg, y0)
    xx=np.linspace(0,H,21)
    err=np.max(np.abs(sol.sol(xx)[0]-u2_exact(xx,H)))
    print(f"H={H}: max|u2_num - u2_exact|={err:.2e}, tol1e-6 -> {'PASS' if err<1e-6 else 'FAIL'}")

print()
print("=== CHECK 5: integration-by-parts identity (quadratic part) ===")
# int (u2'^2 + u2^2/ell^2) dx = [u2 u2']_0^H
for H in [0.5,1.0,2.0]:
    lhs,_=quad(lambda x: u2p_exact(x,H)**2+u2_exact(x,H)**2/ell**2,0,H,limit=200)
    rhs=u2_exact(H,H)*u2p_exact(H,H)-u2_exact(0,H)*u2p_exact(0,H)
    print(f"H={H}: lhs={lhs:.12e} rhs={rhs:.12e} rel={abs(lhs-rhs)/abs(lhs):.2e} -> {'PASS' if abs(lhs-rhs)/abs(lhs)<1e-8 else 'FAIL'}")

print()
print("=== CHECK 6: limits gamma_A, gamma_B ===")
for H in [1e-3,1e-2]:
    beta=H/(2*ell)
    gA=np.tanh(beta); approx_A=beta
    gB=2*np.tanh(beta)-np.tanh(2*beta); approx_B=2*beta**3
    print(f"H={H}: gamma_A={gA:.6e} vs beta={approx_A:.6e}; gamma_B={gB:.6e} vs 2beta^3={approx_B:.6e}")
for H in [10.0,50.0]:
    beta=H/(2*ell)
    print(f"H={H}: gamma_A={np.tanh(beta):.10f} gamma_B={2*np.tanh(beta)-np.tanh(2*beta):.10f} (->1)")

print()
print("=== CHECK 7: G = K^2/E consistency (kI = sqrt(gamma*(1-nu^2))) ===")
for H in [0.5,1.0,2.0]:
    beta=H/(2*ell)
    gA=np.tanh(beta); gB=2*np.tanh(beta)-np.tanh(2*beta)
    kA=np.sqrt((1-nu**2)*gA); kB=np.sqrt((1-nu**2)*gB)
    # reconstruct gamma from kI via Ebar/E = 1/(1-nu^2)
    print(f"H={H}: gamma_A rec={kA**2/(1-nu**2):.12f} orig={gA:.12f}; gamma_B rec={kB**2/(1-nu**2):.12f} orig={gB:.12f}")
