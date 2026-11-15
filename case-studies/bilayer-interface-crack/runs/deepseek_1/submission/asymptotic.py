import numpy as np
from numpy import cos, sin, pi

def F_and_derivs(coef, lam, th):
    A,B,C,D = coef
    a=lam+1.0; b=lam-1.0
    F = A*cos(a*th)+B*sin(a*th)+C*cos(b*th)+D*sin(b*th)
    Fp = -A*a*sin(a*th)+B*a*cos(a*th)-C*b*sin(b*th)+D*b*cos(b*th)
    Fpp = -A*a*a*cos(a*th)-B*a*a*sin(a*th)-C*b*b*cos(b*th)-D*b*b*sin(b*th)
    return F,Fp,Fpp

def stress_disp(coef, lam, th, mu, nu):
    A,B,C,D = coef
    a=lam+1.0; b=lam-1.0
    ca=cos(a*th); sa=sin(a*th); cb=cos(b*th); sb=sin(b*th)
    F = A*ca+B*sa+C*cb+D*sb
    Fp = -A*a*sa+B*a*ca-C*b*sb+D*b*cb
    Fpp = -A*a*a*ca-B*a*a*sa-C*b*b*cb-D*b*b*sb
    Srr = Fpp + a*F
    Stt = lam*a*F
    Srt = -lam*Fp
    kappa = 3.0-4.0*nu
    # analytic displacements
    U = ( -(lam+1.0)*(A*ca+B*sa) + (kappa-lam)*(C*cb+D*sb) )/(2.0*mu)
    Up = ( (lam+1.0)*a*(A*sa-B*ca) + (kappa-lam)*b*(-C*sb+D*cb) )/(2.0*mu)
    V = ( (lam+1.0)*(A*sa-B*ca) + (lam+kappa)*(C*sb-D*cb) )/(2.0*mu)
    Vp = ( (lam+1.0)**2*(A*ca+B*sa) + (lam-1.0)*(lam+kappa)*(C*cb+D*sb) )/(2.0*mu)
    return Srr,Stt,Srt,U,V,Up,Vp

def build_eigen(eps, mu1=3.375, mu2=0.375, nu1=1/3, nu2=1/3):
    lam = 0.5 + 1j*eps
    M = np.zeros((8,8), dtype=complex)
    def rows(th, mu, nu):
        out=[]
        for kind in ['Stt','Srt','ur','ut']:
            r=np.zeros(4,dtype=complex)
            for i in range(4):
                coef=np.zeros(4,dtype=complex); coef[i]=1.0
                Srr,Stt,Srt,U,V,Up,Vp = stress_disp(coef,lam,th,mu,nu)
                val={'Stt':Stt,'Srt':Srt,'ur':U,'ut':V}[kind]
                r[i]=val
            out.append(r)
        return out
    r1=rows(pi,mu1,nu1); r2=rows(-pi,mu2,nu2)
    r0a=rows(0.0,mu1,nu1); r0b=rows(0.0,mu2,nu2)
    M[0,:4]=r1[0]; M[1,:4]=r1[1]
    M[2,4:]=r2[0]; M[3,4:]=r2[1]
    M[4,:4]=r0a[0]; M[4,4:]=-r0b[0]
    M[5,:4]=r0a[1]; M[5,4:]=-r0b[1]
    M[6,:4]=r0a[2]; M[6,4:]=-r0b[2]
    M[7,:4]=r0a[3]; M[7,4:]=-r0b[3]
    U,s,Vh = np.linalg.svd(M)
    v = Vh[-1,:].conj()
    return lam, v

class AuxField:
    """Auxiliary near-tip field with prescribed complex K (task definition)."""
    def __init__(self, eps, mu1=3.375, mu2=0.375, nu1=1/3, nu2=1/3):
        self.lam, self.v = build_eigen(eps, mu1, mu2, nu1, nu2)
        self.mu1, self.mu2, self.nu1, self.nu2 = mu1, mu2, nu1, nu2
        c1=self.v[:4]; c2=self.v[4:]
        Srr,Stt,Srt,U,V,Up,Vp = stress_disp(c1,self.lam,0.0,mu1,nu1)
        self.KF = np.sqrt(2*np.pi)*(Stt+1j*Srt)
        Srr2,Stt2,Srt2,U2,V2,Up2,Vp2 = stress_disp(c2,self.lam,0.0,mu2,nu2)
        self.KF2 = np.sqrt(2*np.pi)*(Stt2+1j*Srt2)
    def _raw(self, x, y, K1, K2):
        K = K1 + 1j*K2
        C = 2.0*K/self.KF
        r = np.hypot(x,y)
        th = np.arctan2(y,x)
        layer = 1 if y >= 0 else 2
        if layer==1:
            coef=self.v[:4]; mu,nu=self.mu1,self.nu1
        else:
            coef=self.v[4:]; mu,nu=self.mu2,self.nu2
        Srr,Stt,Srt,U,V,Up,Vp = stress_disp(coef,self.lam,th,mu,nu)
        return C, r, th, layer, Srr, Stt, Srt, U, V, Up, Vp
    def fields(self, x, y, K1=1.0, K2=0.0):
        C,r,th,layer,Srr,Stt,Srt,U,V,Up,Vp = self._raw(x,y,K1,K2)
        rl = r**self.lam
        rp = r**(self.lam-1.0)
        srr = np.real(C*rp*Srr); stt=np.real(C*rp*Stt); srt=np.real(C*rp*Srt)
        c=cos(th); s=sin(th)
        sxx = srr*c*c + stt*s*s - 2*srt*s*c
        syy = srr*s*s + stt*c*c + 2*srt*s*c
        sxy = (srr-stt)*s*c + srt*(c*c-s*s)
        fx = U*c - V*s; fy = U*s + V*c
        ux = np.real(C*rl*fx); uy=np.real(C*rl*fy)
        fxp = Up*c - U*s - Vp*s - V*c
        fyp = Up*s + U*c + Vp*c - V*s
        uxx = np.real(C*rp*(self.lam*fx*c - fxp*s))
        uxy = np.real(C*rp*(self.lam*fx*s + fxp*c))
        uyx = np.real(C*rp*(self.lam*fy*c - fyp*s))
        uyy = np.real(C*rp*(self.lam*fy*s + fyp*c))
        g = np.array([[uxx,uxy],[uyx,uyy]])
        return np.array([sxx,syy,sxy]), np.array([ux,uy]), g
