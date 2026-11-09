import numpy as np
from scipy.optimize import brentq
import cmath, math

# n=1 (linear incompressible power-law) eigenvalue problem.
# Even solution f = A cos(k1 th) + B cos(k2 th),
#   k^4 - (4-2 s^2) k^2 + s^2 (2-s)^2 = 0
# BCs: f(alpha)=0, f'(alpha)=0  ->  det = 0.
def ks(s):
    b = -(4 - 2*s*s)          # coefficient of k^2 is (4-2s^2) => k^4 + a k^2 + c
    a = (4 - 2*s*s)
    c = s*s*(2-s)**2
    disc = a*a - 4*c
    k1 = cmath.sqrt((a + cmath.sqrt(disc))/2)
    k2 = cmath.sqrt((a - cmath.sqrt(disc))/2)
    return k1, k2

def det(s, alpha):
    k1, k2 = ks(s)
    # f = cos(k1 th)+ r cos(k2 th); f(alpha)=0, f'(alpha)=0
    c1, s1 = cmath.cos(k1*alpha), cmath.sin(k1*alpha)
    c2, s2 = cmath.cos(k2*alpha), cmath.sin(k2*alpha)
    M = np.array([[c1, c2], [-k1*s1, -k2*s2]], dtype=complex)
    return abs(np.linalg.det(M))

if __name__ == "__main__":
    alpha = 2*np.pi/3
    print("n=1 incompressible mode-I eigenvalues s in (0,1), alpha=2pi/3")
    ss = np.linspace(0.01, 0.99, 400)
    vals = [det(s, alpha) for s in ss]
    for i in range(len(ss)-1):
        if vals[i]*vals[i+1] < 0:
            r = brentq(lambda z: det(z, alpha), ss[i], ss[i+1], xtol=1e-13)
            print("  root s =", round(r,6))
    # crack limit alpha->pi
    print("crack limit alpha=pi:")
    for a in [np.pi-1e-4]:
        ss2 = np.linspace(0.01,0.99,400)
        v2=[det(s,a) for s in ss2]
        for i in range(len(ss2)-1):
            if v2[i]*v2[i+1]<0:
                r=brentq(lambda z:det(z,a),ss2[i],ss2[i+1],xtol=1e-13)
                print("  alpha",round(a,4),"root",round(r,6)," (elast/incomp crack 0.5)")