import numpy as np
exec(open('auxfield.py').read().split("if __name__")[0])

def Wrot2(mat,th):
    h=1e-6
    Exx,Eyy,Exy=Eang(mat,th)
    Exx_p=(Eang(mat,th+h)[0]-Eang(mat,th-h)[0])/(2*h)
    Eyy_p=(Eang(mat,th+h)[1]-Eang(mat,th-h)[1])/(2*h)
    Exy_p=(Eang(mat,th+h)[2]-Eang(mat,th-h)[2])/(2*h)
    A=-s*Exy+Exx_p; B=Exy_p+s*Exx; C=Eyy_p+s*Exy; D=-s*Eyy+Exy_p
    co=np.cos(th); si=np.sin(th)
    M=np.array([[co,-si],[si,co]]); rhs=np.array([-(co*A+si*B),-(si*C+co*D)])
    sW,Wp=np.linalg.solve(M,rhs)
    return sW/s
W00={1:Wrot2(1,0.0),2:Wrot2(2,0.0)}
def Wval(mat,th): return Wrot2(mat,th)-W00[mat]
def field_phys(mat,th,r,which='A'):
    Sxx,Syy,Sxy,Exx,Eyy,Exy=complex_field(mat,th,r)
    omega=r**s*Wval(mat,th)
    if which=='A':
        return np.real(Sxx),np.real(Syy),np.real(Sxy),np.real(Exx),np.real(Eyy),np.real(Exy),np.real(omega)
    return np.imag(Sxx),np.imag(Syy),np.imag(Sxy),np.imag(Exx),np.imag(Eyy),np.imag(Exy),np.imag(omega)
def Jcircle(r0=1.0,n=4000,which='A'):
    th=np.linspace(-np.pi,np.pi,n,endpoint=False)+np.pi/n
    vals=np.zeros(n)
    for i,t in enumerate(th):
        mat=1 if t>0 else 2
        sx,sy,sxy,ex,ey,exy,om=field_phys(mat,t,r0,which)
        S=np.array([[sx,sxy],[sxy,sy]]); E=np.array([[ex,exy],[exy,ey]])
        u1=np.array([ex,exy+om]); nvec=np.array([np.cos(t),np.sin(t)])
        W=0.5*np.sum(S*E); tvec=S@nvec
        vals[i]=(W*nvec[0]-np.dot(tvec,u1))*r0
    return np.sum(vals)*2*np.pi/n
for r0 in [0.3,0.6,1.0,2.0,4.0]:
    for which in ['A','B']:
        J=Jcircle(r0,4000,which)
        print("r0",r0,which,"J =",J)
print("expected C|K|^2",Cc*abs(KA)**2)
