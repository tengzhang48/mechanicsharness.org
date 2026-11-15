import numpy as np
exec(open('auxfield.py').read().split("if __name__")[0])

def Wrot2(mat,th):
    h=1e-6
    Exx,Eyy,Exy=Eang(mat,th)
    Exx_p=(Eang(mat,th+h)[0]-Eang(mat,th-h)[0])/(2*h)
    Eyy_p=(Eang(mat,th+h)[1]-Eang(mat,th-h)[1])/(2*h)
    Exy_p=(Eang(mat,th+h)[2]-Eang(mat,th-h)[2])/(2*h)
    A=-s*Exy+Exx_p; B=Exy_p+s*Exx; C=Eyy_p+s*Exy; D=-s*Eyy+Exy_p
    co=np.cos(th); si=np.sin(th); M=np.array([[co,-si],[si,co]])
    rhs=np.array([-(co*A+si*B),-(si*C+co*D)])
    sW,Wp=np.linalg.solve(M,rhs)
    return sW/s
W00={1:Wrot2(1,0.0),2:Wrot2(2,0.0)}
def vals(mat,th,r,which='A'):
    Sxx,Syy,Sxy,Exx,Eyy,Exy=complex_field(mat,th,r)
    om=r**s*(Wrot2(mat,th)-W00[mat])
    if which=='A': return np.real(Exx),np.real(om)
    return np.imag(Exx),np.imag(om)

# integrate from tip along crack faces: theta=pi for material1, theta=-pi for material2
# x=-r, dx=-dr; delta_x(r)=int_0^r (eps1-eps2) dr  (since u_x(-r)-u_x(0)=int_0^{-r} eps dx = int_0^r eps dr)
# similarly delta_y
def jump(r,which='A',n=20000):
    # integrate 0..r with substitution r=u^2 to handle sqrt singularity? use log grid/quad
    # use Gauss-Legendre on [0,r] after substitution r=t^2
    t,w=np.polynomial.legendre.leggauss(n)
    tt=0.5*(t+1)*r; ww=0.5*w*r
    # r = u^2, dr=2u du ; u from 0..sqrt(r)
    u=np.sqrt(tt); dr=2*u*ww
    e1,o1=vals(1,np.pi,tt,which)
    e2,o2=vals(2,-np.pi,tt,which)
    dx=np.sum((e1-e2)*dr); dy=np.sum((o1-o2)*dr)
    return dx,dy

for which in ['A','B']:
    print("field",which)
    for r in [0.01,0.03,0.1,0.3,1.0]:
        dx,dy=jump(r,which)
        # predicted delta = C K r^(1/2+i eps); C = const*(1-2i eps) (real factor)
        Cpred=(1-2j*eps)
        pred=Cpred*(KA if which=='A' else KB)*r**(0.5+1j*eps)
        print(" r",r,"dx,dy",dx,dy," delta_y+i delta_x",dy+1j*dx,"pred(up to real)",pred)
