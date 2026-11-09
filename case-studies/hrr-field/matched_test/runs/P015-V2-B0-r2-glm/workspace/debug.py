import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'knowledge', 'checks'))

E, nu, sy = 70000.0, 0.30, 10.0
Jg, Czz_m = 0.800, -0.005
s2p = np.sqrt(2*np.pi)
T = 4816.47218664297
nue = 0.5 - (1-2*nu)*T/(2*E)
K = np.sqrt(Jg*T)
kappa = (3-nue)/(1+nue)

def f_cart(th):
    c2, s2 = np.cos(th/2), np.sin(th/2)
    s3, c3 = np.sin(3*th/2), np.cos(3*th/2)
    return c2*(1-s2*s3), c2*(1+s2*s3), c2*s2*c3

def u_closed(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    c2, s2 = np.cos(th/2), np.sin(th/2)
    pref = K*(1+nue)/T*np.sqrt(r/s2p)
    return pref*c2*(kappa-1+2*s2**2), pref*s2*(kappa+1-2*c2**2)

rr, aa = 0.05, 0.6
x, y = rr*np.cos(aa), rr*np.sin(aa)
a,b,c = f_cart(aa); sc = K/np.sqrt(2*np.pi*rr)
exx_l = (a-nue*b)*sc/T; eyy_l = (b-nue*a)*sc/T; exy_l = (1+nue)*c*sc/T

for h in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    duxdx = (u_closed(x+h,y)[0]-u_closed(x-h,y)[0])/(2*h)
    print(f"h={h:.0e}: duxdx={duxdx:.6e}  law exx={exx_l:.6e}  diff={duxdx-exx_l:.2e}")

# J decomposition, linear field
def Jparts(R, Wmode, n=200001):
    ang = np.linspace(-np.pi, np.pi, n)
    x, y = R*np.cos(ang), R*np.sin(ang)
    r, th = np.hypot(x,y), np.arctan2(y,x)
    A,B,C = f_cart(th); scl = K/np.sqrt(2*np.pi*r)
    sxx, syy, sxy = A*scl, B*scl, C*scl
    tr = sxx+syy
    We = 0.5*((1+nu)/E*(sxx**2+syy**2+2*sxy**2) - nu/E*tr**2)
    sd = [sxx-tr/3, syy-tr/3, -tr/3, sxy]
    Se = np.sqrt(1.5*(sd[0]**2+sd[1]**2+sd[2]**2+2*sd[3]**2))
    Wp = np.where(Se>sy, 0.5*(1/T-1/E)*(Se**2-sy**2), 0.0)
    if Wmode == 'lin':
        Wpl = 0.5*(1/T-1/E)*Se**2     # linear-hardening plastic part (no sigma_y cutoff)
        W = We + Wpl
    else:
        W = We + Wp
    hR = 1e-6*R
    dux = (u_closed(x+hR,y)[0]-u_closed(x-hR,y)[0])/(2*hR)
    duy = (u_closed(x+hR,y)[1]-u_closed(x-hR,y)[1])/(2*hR)
    tx, ty = sxx*np.cos(ang)+sxy*np.sin(ang), sxy*np.cos(ang)+syy*np.sin(ang)
    JW = np.trapezoid(W*np.cos(ang)*R, ang)
    Jt = np.trapezoid((tx*dux+ty*duy)*R, ang)
    return JW - Jt, JW, Jt

for R in [1e-4, 1e-2, 0.06, 0.2]:
    jlin = Jparts(R, 'lin'); jtrue = Jparts(R, 'true')
    print(f"R={R:7.4f}: J_lin={jlin[0]:.6f} (JW={jlin[1]:.6f}, Jt={jlin[2]:.6f})  "
          f"J_true={jtrue[0]:.6f} (JW={jtrue[1]:.6f}, Jt={jtrue[2]:.6f})")
print("target 0.8; K^2/T =", K**2/T)
