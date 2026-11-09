"""Second route for the proposed T2 output delta(x = 5.00 mm): closed-form profile
u(z) = 1 - (cosh z + cos z)/2 + kappa/2 (cosh z - cos z) + A/2 (sinh z - sin z), with
lambda from 1 + cosh(l)cos(l) - kappa sinh(l) sin(l) = 0 (smallest admissible root) and
u'''(lambda) = 0 for A. Derived here from u'''' = u - 1 with u(0)=u'(0)=0, u''(0)=kappa."""
import math
from scipy.optimize import brentq
E1,h1,E2,h2,sb,dc,G0 = 110e9,0.60e-3,190e9,0.90e-3,0.20e6,0.20e-3,1.40
C = 12/(E1*h1**3) + 12/(E2*h2**3); ell = (dc/(C*sb))**0.25; kap = math.sqrt(2*C*G0)*ell**2/dc
u = lambda z,A: 1-(math.cosh(z)+math.cos(z))/2 + kap/2*(math.cosh(z)-math.cos(z)) + A/2*(math.sinh(z)-math.sin(z))
u3 = lambda z,A: -(math.sinh(z)+math.sin(z))/2 + kap/2*(math.sinh(z)-math.sin(z)) + A/2*(math.cosh(z)+math.cos(z))
Aof = lambda l: ((math.sinh(l)+math.sin(l)) - kap*(math.sinh(l)-math.sin(l)))/(math.cosh(l)+math.cos(l))
f = lambda l: u(l, Aof(l)) - 1
lam = brentq(f, 1.0, 2.0, xtol=1e-15); A = Aof(lam)
print({"lambda": lam, "A": A, "u3_at_L": u3(lam, A), "L_mm": lam*ell*1e3,
       "delta_at_5mm_um": u(5e-3/ell, A)*dc*1e6, "delta_mid_um": u(lam/2, A)*dc*1e6})
