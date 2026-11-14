import numpy as np
from scipy.integrate import quad
from scipy.special import k0 as besselk0

# --- direct real-line evaluation of two reference integrals ---
def F1(X):
    # int_{-inf}^{inf} e^{ipX} (1+ip)^{-1/2} dp, split even/odd
    f = lambda p: (1+1j*p)**-0.5
    re = lambda p: np.real(f(p)); im = lambda p: np.imag(f(p))
    Ic = quad(lambda p: re(p), 0, np.inf, weight='cos', wvar=X, limlst=200)[0]
    Is = quad(lambda p: im(p), 0, np.inf, weight='sin', wvar=X, limlst=200)[0]
    return 2*(Ic - 1j*Is)   # int e^{ipX} f = 2 int_0^inf [re cos - i? ] careful below

# int_{-inf}^{inf} e^{ipX} f(p) dp = 2 int_0^inf [Re f * cos(pX) - Im f * sin(pX)] dp
def F1b(X):
    f = lambda p: (1+1j*p)**-0.5
    Ic = quad(lambda p: np.real(f(p)), 0, np.inf, weight='cos', wvar=X, limlst=200)[0]
    Is = quad(lambda p: np.imag(f(p)), 0, np.inf, weight='sin', wvar=X, limlst=200)[0]
    return 2*(Ic - Is)

def F2(x):
    g = lambda p: 1.0/np.sqrt(p*p+1.0)   # even, real
    return 2*quad(g, 0, np.inf, weight='cos', wvar=x, limlst=200)[0]

print("test1: direct vs exact 2*sqrt(pi)*e^-X X^-1/2")
for X in [0.3, 1.0, 5.0]:
    num = F1b(X); ex = 2*np.sqrt(np.pi)*np.exp(-X)*X**-0.5
    print(f"  X={X}: direct={num:.10f}  exact={ex:.10f}  ratio={num/ex:.10f}")

print("test2: direct vs exact 2*K0(x)")
for x in [0.5, 2.0]:
    num = F2(x); ex = 2*besselk0(x)
    print(f"  x={x}: direct={num:.10f}  exact={ex:.10f}  ratio={num/ex:.10f}")

# --- check cut-side values by near-cut sampling ---
print("\ncut-side sampling, upper cut p = eps + i(1+s), eps=1e-6:")
for s in [0.3, 1.0]:
    eps = 1e-6
    # (1+ip)^{-1/2}
    pR = eps + 1j*(1+s); pL = -eps + 1j*(1+s)
    vR = (1+1j*pR)**-0.5; vL = (1+1j*pL)**-0.5
    print(f"  s={s}: (1+ip)^-1/2  R={vR:.8f}  L={vL:.8f}   formula R={np.exp(-1j*np.pi/2)*s**-0.5:.8f}  L={np.exp(+1j*np.pi/2)*s**-0.5:.8f}")
    # sqrt(p^2+1)
    rR = np.sqrt(pR*pR+1); rL = np.sqrt(pL*pL+1)
    rho = np.sqrt(s*(s+2))
    print(f"        sqrt(p^2+1)  R={rR:.8f}  L={rL:.8f}   formula R={1j*rho:.8f}  L={-1j*rho:.8f}")
print("\nlower cut p = eps - i(1+s):")
for s in [0.3, 1.0]:
    eps = 1e-6
    pR = eps - 1j*(1+s); pL = -eps - 1j*(1+s)
    rR = np.sqrt(pR*pR+1); rL = np.sqrt(pL*pL+1)
    rho = np.sqrt(s*(s+2))
    print(f"  s={s}: sqrt(p^2+1)  R={rR:.8f}  L={rL:.8f}   formula R={-1j*rho:.8f}  L={+1j*rho:.8f}")
    vR = (1+1j*pR)**-0.5; vL = (1+1j*pL)**-0.5
    print(f"        (1+ip)^-1/2  R={vR:.8f}  L={vL:.8f}   formula both={(2+s)**-0.5:.8f}")
