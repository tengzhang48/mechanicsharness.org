import numpy as np
from scipy.integrate import quad

X = 1.0
# f(p) = e^{ipX} (1+ip)^{-1/2} ; exact real-line integral = 2 sqrt(pi) e^{-X} X^{-1/2}
f = lambda p: np.exp(1j*p*X)*(1+1j*p)**-0.5

# real-line integral (validated method)
Ic = quad(lambda p: np.real(np.exp(1j*p*X)*(1+1j*p)**-0.5), 0, np.inf, weight='cos', wvar=X, limlst=300)[0]
Is = quad(lambda p: np.imag(np.exp(1j*p*X)*(1+1j*p)**-0.5), 0, np.inf, weight='sin', wvar=X, limlst=300)[0]
Ireal = 2*(Ic - Is)
print("I real-line =", Ireal, " exact =", 2*np.sqrt(np.pi)*np.exp(-X)*X**-0.5)

# cut integral: sides p = i(1+s) +- eps
eps = 1e-8
g = lambda s, side: np.exp(1j*(side*eps + 1j*(1+s))*X)*(1 + 1j*(side*eps+1j*(1+s)))**-0.5
Icut = quad(lambda s: np.real(g(s,+1)-g(s,-1)), 0, 80, limit=400)[0] \
     + 1j*quad(lambda s: np.imag(g(s,+1)-g(s,-1)), 0, 80, limit=400)[0]
print("I cut (R-L)  =", Icut, " analytic =", -2j*np.pi*np.exp(-X)*X**-0.5)

# big upper arc at radius R, theta in [0.05, pi-0.05] (avoid junction regions)
for R in [50.0, 200.0]:
    th = np.linspace(0.03, np.pi-0.03, 200001)
    arc = np.trapz(f(R*np.exp(1j*th))*1j*R*np.exp(1j*th), th)
    print(f"arc R={R}: {arc}")
# junction pieces near theta=0 and pi, integrated finely
for R in [50.0, 200.0]:
    th0 = np.linspace(0, 0.03, 200001)
    thpi = np.linspace(np.pi-0.03, np.pi, 200001)
    j0 = np.trapz(f(R*np.exp(1j*th0))*1j*R*np.exp(1j*th0), th0)
    jpi = np.trapz(f(R*np.exp(1j*thpi))*1j*R*np.exp(1j*thpi), thpi)
    print(f"junctions R={R}: near0={j0}  nearpi={jpi}  sum={j0+jpi}")
