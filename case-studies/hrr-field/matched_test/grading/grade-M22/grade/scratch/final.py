import numpy as np
from scipy.optimize import brentq

E,nu,sy,J = 70000.0,0.30,10.0,0.800
Czz = 0.005000

nus = lambda T: 0.5+(nu-0.5)*T/E
T = brentq(lambda T: 2*nus(T)*np.sqrt(J/(2*np.pi*T))-Czz, 1e-6, E, xtol=1e-13, rtol=1e-15)
ns = nus(T); K = np.sqrt(J*T)
print("T        = %.6f MPa"%T)
print("nu*      = %.8f"%ns)
print("K        = %.8f N/mm^1.5"%K)
print("K^2/T    = %.10f (=J)"%(K*K/T))

def S_cart(th):
    ch,sh = np.cos(th/2), np.sin(th/2)
    return (ch*(1-sh*np.sin(1.5*th)), ch*(1+sh*np.sin(1.5*th)), ch*sh*np.cos(1.5*th))
def S_pol(th):
    Sxx,Syy,Sxy = S_cart(th); c,s=np.cos(th),np.sin(th)
    return (Sxx*c*c+Syy*s*s+2*Sxy*s*c, Sxx*s*s+Syy*c*c-2*Sxy*s*c, (Syy-Sxx)*s*c+Sxy*(c*c-s*s))
pref = K/np.sqrt(2*np.pi)
th = np.deg2rad(60.0); Srr,Stt,Srt = S_pol(th)
Crr = pref/T*(Srr-ns*Stt); Ctt = pref/T*(Stt-ns*Srr); Cgr = pref/T*2*(1+ns)*Srt
# cross-check by Cartesian strain coefficient transform
Sxx,Syy,Sxy = S_cart(th)
Cxx = pref/T*(Sxx-ns*Syy); Cyy = pref/T*(Syy-ns*Sxx); Cxy = pref/T*(1+ns)*Sxy
c,s=np.cos(th),np.sin(th)
Crr2=Cxx*c*c+Cyy*s*s+2*Cxy*s*c; Ctt2=Cxx*s*s+Cyy*c*c-2*Cxy*s*c; Cgr2=2*((Cyy-Cxx)*s*c+Cxy*(c*c-s*s))
print("theta=60 Crr=%.8f (chk %.8f)"%(Crr,Crr2))
print("theta=60 Ctt=%.8f (chk %.8f)"%(Ctt,Ctt2))
print("theta=60 C_grthet=%.8f (chk %.8f)"%(Cgr,Cgr2))
Cy = K*(3-ns)/(2*T*np.sqrt(np.pi))
print("Cy(90)   = %.8f"%Cy)

# ---- screening interval on theta=0, leading field ----
# max |principal strain| = |Czz(0)|/sqrt(r) = 0.005/sqrt(r)  (C_xx,C_yy smaller)
r_strain = (Czz/0.05)**2                      # >= this
r_see    = K**2/(2*np.pi*(10*sy)**2)          # <= this
r_cap    = 0.20
print("strain screen  r >= %.8f mm"%r_strain)
print("sigma_e screen r <= %.8f mm"%r_see)
print("r<=0.2 cap active?", r_cap<r_see)
print("INTERVAL = [%.6f, %.6f] mm"%(r_strain, min(r_see,r_cap)))
