import numpy as np
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])

# extract potentials from Airy coefficients: f = a e^{is th}+ b e^{-is th}+ c e^{i(s+2)th}+ d e^{-i(s+2)th}
a,b,c,d = c1
A = a + np.conj(b)
C = c + np.conj(d)
lam = s+1
B = (s+2)*C
print("A=",A," B=",B," lam=",lam)

def pot_field(cAiry, x, y, mat=1):
    # use potentials
    A = cAiry[0] + np.conj(cAiry[1])
    C = cAiry[2] + np.conj(cAiry[3])
    B = (s+2)*C
    z = x+1j*y
    zc = np.conj(z)
    phi = A*z**lam
    psi = B*z**lam
    Phi = lam*A*z**(lam-1)
    Phip = lam*(lam-1)*A*z**(lam-2)
    Phic = np.conj(lam)*np.conj(A)*zc**(np.conj(lam)-1)
    Phipc = np.conj(lam)*(np.conj(lam)-1)*np.conj(A)*zc**(np.conj(lam)-2)
    psic = np.conj(B)*zc**np.conj(lam)
    Psi = lam*B*z**(lam-1)
    F = Phi + Phic + z*Phipc + np.conj(Psi)   # sig_yy - i sig_xy
    tr = 2*(Phi+Phic)
    syy = np.real(F); sxy = -np.imag(F)
    sxx = np.real(tr) - syy
    E,nup,mu = (Estar1,nuprime,mu1) if mat==1 else (Estar2,nuprime,mu2)
    D = (kappa*phi - z*np.conj(Phi) - np.conj(psi))/(2*mu)
    return sxx,syy,sxy,np.real(D),np.imag(D)

# compare stresses with Airy field
def airy_stress(cAiry,x,y):
    r=np.hypot(x,y); th=np.arctan2(y,x)
    f=cAiry@basis(th,s); fp=cAiry@dbasis(th,s); fpp=cAiry@ddbasis(th,s)
    Prr=(s+2)*f+fpp; Ptt=(s+2)*(s+1)*f; Prt=-(s+1)*fp
    rs=r**s
    Sxx=rs*(Prr*np.cos(th)**2+Ptt*np.sin(th)**2-2*Prt*np.sin(th)*np.cos(th))
    Syy=rs*(Prr*np.sin(th)**2+Ptt*np.cos(th)**2+2*Prt*np.sin(th)*np.cos(th))
    Sxy=rs*((Prr-Ptt)*np.sin(th)*np.cos(th)+Prt*(np.cos(th)**2-np.sin(th)**2))
    return np.real(Sxx),np.real(Syy),np.real(Sxy)

for mat,cA in [(1,c1),(2,c2)]:
    print("material",mat)
    for (x,y) in [(0.5,0.3),(1.0,0.2),(-0.7,0.4),(0.2,-0.6),(-1.0,-0.5),(2.0,1.0)]:
        pa=pot_field(cA,x,y,mat)
        aa=airy_stress(cA,x,y)
        print("  point",x,y)
        print("    pot ",np.round(pa[:3],8))
        print("    airy",np.round(aa,8))
