import numpy as np
exec(open('auxfield.py').read().split("if __name__")[0])

def Wval(mat,th):
    sW,Wp,res=Wrot(mat,th,1.0)
    return sW

def field_phys(mat,th,r,which='A'):
    Sxx,Syy,Sxy,Exx,Eyy,Exy=complex_field(mat,th,r)
    W=Wval(mat,th)*r**s if False else None
    # Wrot returns W at r=1; omega = r^s W(theta)
    W0=Wval(mat,th)
    omega=r**s*W0
    if which=='A':
        sx,sy,sxy=np.real(Sxx),np.real(Syy),np.real(Sxy)
        ex,ey,exy=np.real(Exx),np.real(Eyy),np.real(Exy)
        om=np.real(omega)
    else:
        sx,sy,sxy=np.imag(Sxx),np.imag(Syy),np.imag(Sxy)
        ex,ey,exy=np.imag(Exx),np.imag(Eyy),np.imag(Exy)
        om=np.imag(omega)
    return sx,sy,sxy,ex,ey,exy,om

def Jcircle(r0=1.0,n=200000,which='A'):
    th=np.linspace(-np.pi,np.pi,n,endpoint=False)+np.pi/n
    # midpoint, not hitting crack faces
    integ=0.0
    vals=[]
    for t in th:
        mat=1 if t>0 else 2
        sx,sy,sxy,ex,ey,exy,om=field_phys(mat,t,r0,which)
        # Cartesian stress tensor, strains, u_i1
        S=np.array([[sx,sxy],[sxy,sy]])
        E=np.array([[ex,exy],[exy,ey]])
        u1=np.array([ex, exy+om])  # u_{x,1}, u_{y,1}
        n=np.array([np.cos(t),np.sin(t)])
        W=0.5*np.sum(S*E)
        tvec=S@n
        val=(W*n[0]-np.dot(tvec,u1))*r0
        integ+=val
    return integ*2*np.pi/n

for which in ['A','B']:
    J=Jcircle(1.0,200000,which)
    print(which,"J =",J)
print("C|KA|^2",Cc*abs(KA)**2)
