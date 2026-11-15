import numpy as np
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])

# complex Airy field, material 1 or 2, returns complex Cartesian stresses and strains (not real part)
def complex_field(mat, th, r=1.0):
    c = c1 if mat==1 else c2
    E,nup,mu = (Estar1,nuprime,mu1) if mat==1 else (Estar2,nuprime,mu2)
    f=c@basis(th,s); fp=c@dbasis(th,s); fpp=c@ddbasis(th,s)
    Prr=(s+2)*f+fpp; Ptt=(s+2)*(s+1)*f; Prt=-(s+1)*fp
    rs=r**s
    Sxx=rs*(Prr*np.cos(th)**2+Ptt*np.sin(th)**2-2*Prt*np.sin(th)*np.cos(th))
    Syy=rs*(Prr*np.sin(th)**2+Ptt*np.cos(th)**2+2*Prt*np.sin(th)*np.cos(th))
    Sxy=rs*((Prr-Ptt)*np.sin(th)*np.cos(th)+Prt*(np.cos(th)**2-np.sin(th)**2))
    Exx=(Sxx-nup*Syy)/E
    Eyy=(Syy-nup*Sxx)/E
    Exy=Sxy/(2*mu)
    return Sxx,Syy,Sxy,Exx,Eyy,Exy

# angular functions E_ij(theta) = r^{-s} eps_ij at r=1
def Eang(mat,th):
    return complex_field(mat,th,1.0)[3:]

def Emat(mat,th,r):
    Exx,Eyy,Exy=Eang(mat,th)
    return r**s*np.array([[Exx,Exy],[Exy,Eyy]],dtype=complex)

def Wrot(mat,th,r):
    # solve for omega = r^s W(theta) from compatibility
    h=1e-6
    Exx,Eyy,Exy=Eang(mat,th)
    Exx_p=(Eang(mat,th+h)[0]-Eang(mat,th-h)[0])/(2*h)
    Eyy_p=(Eang(mat,th+h)[1]-Eang(mat,th-h)[1])/(2*h)
    Exy_p=(Eang(mat,th+h)[2]-Eang(mat,th-h)[2])/(2*h)
    A=-s*Exy+Exx_p
    B=Exy_p+s*Exx
    C=-s*Eyy+Exy_p
    D=-Eyy_p+s*Exy
    # equations: cos*(sW+A) + sin*(-Wp+B)=0 ; sin*(sW+C)+cos*(Wp+D)=0
    # M [sW, Wp]^T = rhs
    co=np.cos(th); si=np.sin(th)
    M=np.array([[co,-si],[si,co]])
    rhs=np.array([-(co*A+si*B), -(si*C+co*D)])
    sW,Wp=np.linalg.solve(M,rhs)
    return sW/s, Wp, (M@np.array([sW,Wp])-rhs)

if __name__=='__main__':
    for mat in [1,2]:
        print("material",mat)
        for th in [0.2,0.7,1.5,2.5,-0.5,-1.5,-2.5]:
            sW,Wp,res=Wrot(mat,th,1.0)
            print(" th",th,"W",sW,"Wp",Wp,"res",res)
