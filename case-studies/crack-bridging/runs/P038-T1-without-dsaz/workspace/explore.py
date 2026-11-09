import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# ---- Parameters ----
E = 130e9
nu = 0.28
Ep = E/(1-nu**2)          # plane strain modulus
a = 2.40e-3               # original slit half length (m)
t = 160e6                 # connector closing traction (Pa)
dc = 9.00e-6              # critical full gap (m)
G0 = 22.0                 # intrinsic front work (J/m^2)

# ---- Exact K at outer front, central crack half-length c ----
# p(x)=Sigma on [0,a], Sigma-t on [a,c]
def K_closed(c, S):
    return np.sqrt(c/np.pi)*(np.pi*S - t*(np.pi - 2*np.arcsin(a/c)))

# numeric verification of K formula
def K_num(c, S):
    f = lambda x: (S if x < a else S-t)/np.sqrt(1-(x/c)**2)
    val,_ = quad(f, 0, c, points=[a], limit=400)
    return (2/np.sqrt(np.pi*c))*val

# ---- Exact full gap at crack mouth x=a ----
def delta_a_closed(c, S):
    return (4/Ep)*(S*np.sqrt(c**2-a**2) - t*(c-a))

# numeric verification of general delta(x) formula at x=a
def delta_num(x, c, S):
    def inner(s):
        lo = max(x, a)
        if lo >= c: return 0.0
        val,_ = quad(lambda ss: np.sqrt(ss**2-a**2)/np.sqrt(ss**2-x**2), lo, c,
                     points=[x] if lo==x and x>a else None, limit=400)
        return val
    # delta = (4/Ep)( S*sqrt(c^2-x^2) - t*inner )
    return (4/Ep)*(S*np.sqrt(c**2-x**2) - t*inner(x))

if __name__ == "__main__":
    # checks
    c = a + 1e-3
    S = 100e6
    print("K closed vs num:", K_closed(c,S), K_num(c,S))
    print("delta_a closed vs num:", delta_a_closed(c,S), delta_num(a,c,S))
    # L=0 check: K should be Sigma*sqrt(pi*a)
    print("L=0 K check:", K_closed(a, 50e6), 50e6*np.sqrt(np.pi*a))