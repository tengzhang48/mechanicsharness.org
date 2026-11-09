"""Validate corrected RHS: at n=1 the Williams Mode I Airy eigenfunction must
be a solution (residual ~0). Also probe shooting behaviour for n=6.5."""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp

# ---- n=1 check ----
def williams(t):
    # Airy eigenfunction for Mode I lambda=1/2 with our k=3/2:
    # stress ~ r^{-1/2}. Standard Irwin field angular functions (unit normaliz):
    # sigtilde_rr = cos(t/2)(1 - sin(t/2)sin(3t/2)) etc.  Build phi from Airy:
    # A = k phi + phi'' = sigtilde_rr ; B = k(k-1)phi = sigtilde_tt ; phi = B/(k(k-1)).
    k = 1.5
    stt = np.cos(t/2)*(1 + np.sin(t/2)*np.sin(3*t/2))
    # B = k(k-1) phi => phi = stt/(k(k-1)) = stt/0.75
    return stt/0.75

def williams_d(t, m):
    k = 1.5
    h = 1e-4
    if m == 0:
        return williams(t)
    if m == 1:
        return (williams(t+h)-williams(t-h))/(2*h)
    if m == 2:
        return (williams(t+h)-2*williams(t)+williams(t-h))/h**2
    if m == 3:
        return (williams(t+2*h)-2*williams(t+h)+2*williams(t-h)-williams(t-2*h))/(2*h**3)
    if m == 4:
        return (williams(t+2*h)-4*williams(t+h)+6*williams(t)-4*williams(t-h)+williams(t-2*h))/h**4

print("=== n=1 ridge check: residual phi'''' - rhs_phi'''' ===")
for t in [0.3, 1.0, 2.0, 2.9]:
    Y = [williams_d(t,0), williams_d(t,1), williams_d(t,2), williams_d(t,3)]
    r = H.rhs(t, np.array(Y), n=1.0, p=0.5, k=1.5)
    print(f"  t={t}: phi''''_fd={williams_d(t,4):.5f}  rhs={r[3]:.5f}  diff={williams_d(t,4)-r[3]:.2e}")