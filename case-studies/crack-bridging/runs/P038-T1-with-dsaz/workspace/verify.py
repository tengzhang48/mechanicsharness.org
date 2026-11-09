import numpy as np
from scipy.optimize import brentq
E=130e9; nu=0.28; Ep=E/(1-nu**2); s0=160e6; dc=9e-6; G0=22.0; a=2.4e-3
K0=np.sqrt(Ep*G0); ell=Ep*dc/s0

# 1) Consistency: eq(46) with K_tip=0, delta_t=dc must reproduce eq(47).
def sbar47(eta): return s0*(2/np.pi)*np.arccos(1/eta)
def dt47(eta):
    r = sbar47(eta)/s0 - (2/np.pi)*np.arccos(1/eta)  # = 0
    return s0*a/Ep*( (8/np.pi)*np.log(eta) )
# solve dt47=dc -> eq(47) predicts
eta47 = brentq(lambda e: dt47(e)-dc, 1.001, 5, xtol=1e-14)
S47 = sbar47(eta47)/s0
print(f'Gamma0=0 limit: eta={eta47:.9f}  S/s0={S47:.9f}')
# eq(47): (8/pi)*alpha*ln(sec(pi/2 * S/s0)) = 1
alpha=a/ell
from math import pi, acos, exp, cos, log
S47b = 2/pi*acos(exp(-pi/(8*alpha)))
print(f'  eq(47) closed form S/s0={S47b:.9f}  (match={abs(S47-S47b):.2e})')

# 2) Recompute the answer with alternative (sigma_bar, eta) both unknown via 2 eqs
from scipy.optimize import fsolve
def eqs(v):
    S, eta = v
    r = S/s0 - (2/np.pi)*np.arccos(1/eta)
    e45 = s0*a/Ep*(4*np.sqrt(eta**2-1)*r + (8/np.pi)*np.log(eta)) - dc
    e46 = s0*np.sqrt(np.pi*a)*np.sqrt(eta)*r - K0
    return [e45, e46]
sol = fsolve(eqs, [120e6, 2.7], full_output=True)
S, eta = sol[0]
print(f'2-eq solve: Sigma={S/1e6:.6f} MPa  eta={eta:.9f}  L={(eta-1)*a*1e3:.6f} mm')
print(f'  residual={sol[1]["fvec"]}')
