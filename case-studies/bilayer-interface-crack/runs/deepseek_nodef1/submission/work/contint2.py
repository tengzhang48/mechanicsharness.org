import numpy as np
exec(open('auxfield.py').read().split("if __name__")[0])

def field_phys(mat,th,r,which='A'):
    Sxx,Syy,Sxy,Exx,Eyy,Exy=complex_field(mat,th,r)
    W0=Wrot(mat,th,1.0)[0]
    omega=r**s*W0
    if which=='A':
        return np.real(Sxx),np.real(Syy),np.real(Sxy),np.real(Exx),np.real(Eyy),np.real(Exy),np.real(omega)
    else:
        return np.imag(Sxx),np.imag(Syy),np.imag(Sxy),np.imag(Exx),np.imag(Eyy),np.imag(Exy),np.imag(omega)

def Jcircle(r0=1.0,n=4000,which='A'):
    th=np.linspace(-np.pi,np.pi,n,endpoint=False)+np.pi/n
    vals=np.zeros(n)
    for i,t in enumerate(th):
        mat=1 if t>0 else 2
        sx,sy,sxy,ex,ey,exy,om=field_phys(mat,t,r0,which)
        S=np.array([[sx,sxy],[sxy,sy]]); E=np.array([[ex,exy],[exy,ey]])
        u1=np.array([ex, exy+om]); nvec=np.array([np.cos(t),np.sin(t)])
        W=0.5*np.sum(S*E); tvec=S@nvec
        vals[i]=(W*nvec[0]-np.dot(tvec,u1))*r0
    return np.sum(vals)*2*np.pi/n
for r0 in [0.3,0.6,1.0,2.0]:
    for which in ['A','B']:
        J=Jcircle(r0,4000,which)
        print("r0",r0,which,"J =",J, " J/C",J/Cc)
print("|KA|^2",abs(KA)**2,"|KB|^2",abs(KB)**2)
