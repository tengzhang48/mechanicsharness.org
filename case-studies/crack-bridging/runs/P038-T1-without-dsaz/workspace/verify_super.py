import numpy as np
from scipy.integrate import quad
import explore

Ep = explore.Ep; a = explore.a; t = explore.t

# Independent SUPERPOSITION check of K and delta(a).
# Known closed forms for a central crack of half-length c in an infinite plate:
#   (i)  uniform pressure p0 on the WHOLE crack [0,c]:  K = p0*sqrt(pi*c),
#        delta(x) = (4 p0/Ep)*sqrt(c^2 - x^2).
#   (ii) uniform pressure p0 on a PART [0,b] (b<c):     K = 2 p0 sqrt(c/pi) asin(b/c),
#        delta(x) (x<=b) = (4 p0/Ep)*( sqrt(c^2-x^2) - ... ) ; at x=b:
#        delta(b) = (4 p0/Ep)*( b - c + sqrt(c^2-b^2) )   [derived via CRC integral].
# Our load: p = Sigma on [0,a], Sigma - t on [a,c].
#   = (Sigma - t) uniform on [0,c]  +  t uniform on [0,a].
# Therefore:
#   K(c,Sigma) = (Sigma - t)*sqrt(pi*c) + 2 t sqrt(c/pi) asin(a/c)
#   delta(a)   = (4/Ep)[ (Sigma-t)*sqrt(c^2-a^2) + t*(a - c + sqrt(c^2-a^2)) ]

def K_super(c, S):
    return (S - t)*np.sqrt(np.pi*c) + 2*t*np.sqrt(c/np.pi)*np.arcsin(a/c)

def K_closed(c, S):   # from explore (independently validated by numeric quadrature)
    return explore.K_closed(c, S)

def delta_part_b_numeric(b, c, p0):
    # direct CRC integral for uniform p0 on [0,b], evaluated at x=b
    def inner(s):
        if s <= b:
            q = p0*s
        else:
            q = p0*(s - np.sqrt(s*s - b*b))
        return q/np.sqrt(s*s - b*b)
    val,_ = quad(inner, b, c, limit=400, points=[b+1e-12])
    return (4.0/Ep)*val

if __name__=="__main__":
    c = a + 34.324304e-3
    S = 158.524900e6
    print("K superposition = %.10e" % K_super(c,S))
    print("K closed        = %.10e" % K_closed(c,S))
    # delta(a) superposition
    d1 = (4.0/Ep)*((S-t)*np.sqrt(c*c-a*a))
    d2 = delta_part_b_numeric(a, c, t)
    print("delta(a) super = %.10e m" % (d1+d2))
    print("delta(a) closed= %.10e m" % explore.delta_a_closed(c,S))
    print("delta(a) target= %.10e" % explore.dc)

    # also scan: at the solved L*, does G=G0 and delta=dc hold via superposition?
    print("G from super = %.10e J/m2" % (K_super(c,S)**2/Ep))