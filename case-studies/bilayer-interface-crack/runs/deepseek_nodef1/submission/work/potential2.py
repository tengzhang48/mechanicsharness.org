import numpy as np
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])

def potentials_from_airy(cAiry):
    a,b,c,d = cAiry
    A = 2*a
    B = 2*np.conj(b) - 2*np.conj(a)
    C = 2*c
    D = 2*np.conj(d) - 2*np.conj(c)
    return A,B,C,D

def pot_stress(cAiry,x,y):
    A,B,C,D = potentials_from_airy(cAiry)
    z=x+1j*y; zc=np.conj(z)
    lam=s+1
    # phi = A z^lam + B z^conj(lam); chi = C z^(lam+1)+D z^(conj(lam)+1)
    phi = A*z**lam + B*z**np.conj(lam)
    Phi = lam*A*z**(lam-1) + np.conj(lam)*B*z**(np.conj(lam)-1)
    Phip = lam*(lam-1)*A*z**(lam-2) + np.conj(lam)*(np.conj(lam)-1)*B*z**(np.conj(lam)-2)
    chi = C*z**(lam+1)+D*z**(np.conj(lam)+1)
    psi = (lam+1)*C*z**lam + (np.conj(lam)+1)*D*z**np.conj(lam)
    Psi = (lam+1)*lam*C*z**(lam-1) + (np.conj(lam)+1)*np.conj(lam)*D*z**(np.conj(lam)-1)
    # conjugates: bar Phi = conj(A)*conj(lam)*zc^(conj(lam)-1)+...
    Phic = np.conj(lam)*np.conj(A)*zc**(np.conj(lam)-1) + lam*np.conj(B)*zc**(lam-1)
    Phipc = np.conj(lam)*(np.conj(lam)-1)*np.conj(A)*zc**(np.conj(lam)-2) + lam*(lam-1)*np.conj(B)*zc**(lam-2)
    Psic = np.conj(lam+1)*np.conj(C)*zc**np.conj(lam) + (lam+1)*np.conj(D)*zc**lam
    F = Phi + Phic + z*Phipc + np.conj(Psi)  # sig_yy - i sig_xy
    tr = 2*(Phi+Phic)
    syy=np.real(F); sxy=-np.imag(F); sxx=np.real(tr)-syy
    return sxx,syy,sxy,phi,psi,Phi,Phic

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
        ps=pot_stress(cA,x,y)[:3]; aa=airy_stress(cA,x,y)
        print("  ",x,y,"pot",np.round(ps,6),"airy",np.round(aa,6),"diff",np.round(np.array(ps)-np.array(aa),8))
