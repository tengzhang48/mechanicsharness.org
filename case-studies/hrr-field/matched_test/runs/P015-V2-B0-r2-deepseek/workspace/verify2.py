import numpy as np
from scipy.optimize import brentq

E=70000.0; nu=0.30; sy=10.0; J=0.8; Czz=-0.005
P_of = lambda T: 2*nu/E + 1/T - 1/E
Czz_of = lambda T: -P_of(T)*np.sqrt(J*T)/np.sqrt(2*np.pi)
T = brentq(lambda T: Czz_of(T)-Czz, 1, E-1e-6)
K = np.sqrt(J*T)
D=(1+nu)/E+1.5*(1/T-1/E); B=(1-2*nu)/E
nu_eff=(D-B)/(2*D+B); E_eff=T
G_eff=E_eff/(2*(1+nu_eff)); kappa=(3-nu_eff)/(1+nu_eff)

def u_cart(x,y):
    r=np.hypot(x,y); th=np.arctan2(y,x); c=np.cos(th/2); s=np.sin(th/2)
    fac=K/(2*G_eff)*np.sqrt(r/(2*np.pi))
    return fac*c*(kappa-1+2*s*s), fac*s*(kappa+1-2*c*c)

# numeric strain at theta=60 via central differences in Cartesian, then rotate to polar
th=np.radians(60.0); r0=1e-3
x0,y0=r0*np.cos(th),r0*np.sin(th)
h=1e-7
uxp,uyp=u_cart(x0+h,y0); uxm,uym=u_cart(x0-h,y0)
duydx=(uyp-uym)/(2*h)
uxp,uyp=u_cart(x0,y0+h); uxm,uym=u_cart(x0,y0-h)
duydy=(uyp-uym)/(2*h)
duxdx=(u_cart(x0+h,y0)[0]-u_cart(x0-h,y0)[0])/(2*h)
duxdy=(u_cart(x0,y0+h)[0]-u_cart(x0,y0-h)[0])/(2*h)
# rotate strain tensor to polar
c,s=np.cos(th),np.sin(th)
exx,eyy,exy=duxdx,duydy,0.5*(duxdy+duydx)
err=exx*c*c+eyy*s*s+2*exy*s*c
ett=exx*s*s+eyy*c*c-2*exy*s*c
ert=-(exx-eyy)*s*c+exy*(c*c-s*s)
print("numeric at theta=60: C_rr=%.9f C_tt=%.9f C_grth=%.9f" %
      (np.sqrt(r0)*err, np.sqrt(r0)*ett, np.sqrt(r0)*2*ert))

# analytic
c,s=np.cos(th/2),np.sin(th/2)
f_rr=c*(1+s*s); f_tt=c*(1-s*s); f_rt=s*c*c
A=K/(E_eff*np.sqrt(2*np.pi))
print("analytic             C_rr=%.9f C_tt=%.9f C_grth=%.9f" %
      (A*(f_rr-nu_eff*f_tt), A*(f_tt-nu_eff*f_rr), 2*(1+nu_eff)*A*f_rt))

# C_y(90) numeric
ux,uy=u_cart(0.0,1.0)
print("C_y(90) closed=%.9f  num=%.9f" % (uy, u_cart(0,1e-6)[1]/np.sqrt(1e-6)))

# independent verification of effective-modulus reduction: full law at huge stress vs effective Hooke
def full_strain(sig):
    # uniaxial-like plane-stress state sigma_xx=sigma_yy=sig, sigma_zz=0, sigma_e=sig
    sxx=sig; syy=sig; szz=0.0; skk=sxx+syy+szz
    # deviator components
    m=skk/3
    dxx,dzz=sxx-m, szz-m
    se=sig
    active = se>sy
    fac=1.5*(1/T-1/E)*(1-sy/se) if active else 0.0
    exx=(1+nu)/E*sxx - nu/E*skk + fac*dxx
    ezz=(1+nu)/E*szz - nu/E*skk + fac*dzz
    return exx,ezz
sig=1e7
exx,ezz=full_strain(sig)
print("full law at sig=1e7: eps_xx=%.12e eps_zz=%.12e ratio ezz/exx=%.10f"%(exx,ezz,ezz/exx))
print("effective:          eps_xx=%.12e eps_zz=%.12e ratio=%.10f"%((1-nu_eff)*sig/E_eff,-2*nu_eff*sig/E_eff,(-2*nu_eff)/(1-nu_eff)))
