import numpy as np
E, nu, sy = 70000.0, 0.30, 10.0
Jg = 0.800
s2p = np.sqrt(2*np.pi)
T = 4816.47218664297
nue = 0.5 - (1-2*nu)*T/(2*E)
K = np.sqrt(Jg*T)
kappa = (3-nue)/(1+nue)
print(f"T={T:.6f} nue={nue:.8f} K={K:.8f} kappa={kappa:.8f} fac={K/(s2p*T):.8e} A'={K*(1+nue)/(s2p*T):.8e}")

def f_cart(th):
    c2, s2 = np.cos(th/2), np.sin(th/2)
    s3, c3 = np.sin(3*th/2), np.cos(3*th/2)
    return c2*(1-s2*s3), c2*(1+s2*s3), c2*s2*c3

rr, aa = 0.05, 0.6
x, y = rr*np.cos(aa), rr*np.sin(aa)
a,b,c = f_cart(aa); sc = K/np.sqrt(2*np.pi*rr)
exx_l = (a-nue*b)*sc/T
print(f"law exx = {exx_l:.8e}")

def u_closed(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    c2, s2 = np.cos(th/2), np.sin(th/2)
    pref = K*(1+nue)/T*np.sqrt(r/s2p)
    return pref*c2*(kappa-1+2*s2**2), pref*s2*(kappa+1-2*c2**2)

h = 1e-6
up, um = u_closed(x+h, y), u_closed(x-h, y)
print(f"ux(x+h,y)={up[0]:.12e}  ux(x-h,y)={um[0]:.12e}")
print(f"FD duxdx = {(up[0]-um[0])/(2*h):.8e}")

# analytic partials of u_closed at this point
r0, th0 = rr, aa
c2, s2 = np.cos(th0/2), np.sin(th0/2)
G = kappa-1+2*s2**2
pref = K*(1+nue)/T*np.sqrt(r0/s2p)
dpx_dr = pref/(2*r0)
dGdth = -s2/2*G + c2*(2*s2*c2)
dup_dth = pref*dGdth
dup_dr = dpx_dr*c2*G
cx, sy = np.cos(th0), np.sin(th0)
print(f"analytic duxdx = {dup_dr*cx - dup_dth*sy/r0:.8e}")
print(f"parts: dup_dr*cx = {dup_dr*cx:.6e}, dup_dth*sy/r0 = {dup_dth*sy/r0:.6e}")

# Now the integrated displacement (from strains) and its FD
N = 200001
tg = np.linspace(0.0, np.pi, N)
Frr, Ftt, Frt = [], [], []
for th in tg:
    A,B,C = f_cart(th); co, si = np.cos(th), np.sin(th)
    Frr.append(A*co**2 + B*si**2 + 2*C*si*co)
    Ftt.append(A*si**2 + B*co**2 - 2*C*si*co)
    Frt.append((B-A)*si*co + C*(co**2-si**2))
Frr, Ftt, Frt = map(np.array, (Frr, Ftt, Frt))
fac = K/(s2p*T)
err = (Frr - nue*Ftt)*fac
ett = (Ftt - nue*Frr)*fac
ert = (1+nue)*Frt*fac
Ur = 2*err
dUt = ett - Ur
Ut = np.concatenate([[0.0], np.cumsum(0.5*(dUt[1:]+dUt[:-1])*np.diff(tg))])
Ur_s = lambda th: np.interp(th, tg, Ur)
Ut_s = lambda th: np.interp(th, tg, Ut)
def u_int(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    th = np.where(th < 0, -th, th)          # mode I: u_x even, u_y odd
    ur, ut = np.sqrt(r)*Ur_s(th), np.sqrt(r)*Ut_s(th)
    sgn = np.where(np.arctan2(y, x) < 0, -1.0, 1.0)
    return ur*np.cos(th) - sgn*ut*np.sin(th), sgn*(ur*np.sin(th) + ut*np.cos(th))
up2, um2 = u_int(x+h, y), u_int(x-h, y)
print(f"integrated field FD duxdx = {(up2[0]-um2[0])/(2*h):.8e}")
up3, um3 = u_int(x, y+h), u_int(x, y-h)
print(f"integrated field FD duydy = {(up3[1]-um3[1])/(2*h):.8e}  law eyy = {(b-nue*a)*sc/T:.8e}")
# compare closed form vs integrated polar components at several angles
for th in [0.3, 0.6, 1.2, np.pi/2, 2.5]:
    urc = K*(1+nue)/T/np.sqrt(s2p)*np.cos(th/2)*(kappa-1+2*np.sin(th/2)**2)
    utc = -K*(1+nue)/T/np.sqrt(s2p)*np.sin(th/2)*(kappa+1-2*np.cos(th/2)**2)
    print(f"th={th:.3f}: Ur_int={Ur_s(th):.6e} Ur_cf={urc:.6e} | Ut_int={Ut_s(th):.6e} Ut_cf={utc:.6e}")
