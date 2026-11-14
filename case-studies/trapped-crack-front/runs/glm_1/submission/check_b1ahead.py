"""Validate b1_ahead closed forms against direct quadrature of the lower-cut formula:
T1^{ij}(|x|) = -(i/2pi) int_0^inf (fR-fL) e^{-(1+s)|x|} ds,
fR-fL = -2[(i/rho)Phat_ij + i rho Qhat_ij] * sqrt(pi) (2+s)^{-1/2},
rho = sqrt(s(s+2)); Phat/Qhat per the M-entry decompositions (checked vs M_closed).
"""
import numpy as np
from scipy.integrate import quad
from m_stiffness import M_closed
from first_order import b1_ahead, CC, NU

def Phat_Qhat(entry, s):
    rho = np.sqrt(s*(s+2.0))
    if entry == 'xx':
        return CC*NU*np.ones_like(s), CC*np.ones_like(s)
    if entry == 'zz':
        return -CC*NU*np.ones_like(s), CC*(1-NU)*np.ones_like(s)
    if entry == 'xz':
        return 1j*CC*NU*(1.0+s), np.zeros_like(s)
    raise KeyError

def T1_cut_quad(entry, x):
    def g(s):
        P, Q = Phat_Qhat(entry, np.array([s]))
        rho = np.sqrt(s*(s+2.0))
        fRmL = -2*((1j/rho)*P[0] + 1j*rho*Q[0])*np.sqrt(np.pi)*(2.0+s)**-0.5
        return fRmL*np.exp(-(1.0+s)*x)
    val = quad(lambda s: np.real(g(s)), 0, 300, limit=500)[0] \
        + 1j*quad(lambda s: np.imag(g(s)), 0, 300, limit=500)[0]
    return -1j/(2*np.pi)*val

# cross-check decomposition Phat/Qhat against the M formula on the cut sides
def Mentry(entry, p, r):
    cc = CC
    if entry == 'xx': return cc*((1-NU)*r + NU*p*p/r)
    if entry == 'zz': return cc*((1-NU)*r + NU/r)
    if entry == 'xz': return cc*(NU*p/r)
for entry in ['xx', 'zz', 'xz']:
    s = 0.7; rho = np.sqrt(s*(s+2)); p = -1j*(1+s)
    P, Q = Phat_Qhat(entry, np.array([s]))
    Mminus = Mentry(entry, p, +1j*rho)   # r = +i rho side
    Mdec = (1j/rho)*P[0] + 1j*rho*Q[0]
    print(f"decomp check {entry}: M(r=+i rho)={Mminus:.10f}  decomp={Mdec:.10f}")

print("\nahead closed forms vs cut quadrature:")
for entry in ['xx', 'zz', 'xz']:
    for x in [0.01, 0.1, 0.5, 2.0]:
        v1 = b1_ahead(entry, x); v2 = T1_cut_quad(entry, x)
        print(f"  {entry} x={x:5.2f}: closed={v1:+.8f}  quad={v2:+.8f}  diff={v1-v2:+.1e}")
