import numpy as np, fem_j as FJ
ell=1.0
u,xs,ys,nxn,nyn,nid,shape,gauss=FJ.solve(ell,400,80,20.0,10.0)
U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
dx=xs[1]-xs[0]; dy=ys[1]-ys[0]
Ds=1.0/(1-0.3**2)
def stress(i,j):
    i1=min(i+1,nxn-1); i0=max(i-1,0); j1=min(j+1,nyn-1); j0=max(j-1,0)
    u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0]); u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
    u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0]); u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
    e11=u1x;e22=u2y;e12=0.5*(u1y+u2x);tr=e11+e22
    s11=1.0+Ds*((1-0.3)*e11+0.3*tr); s22=1.0+Ds*((1-0.3)*e22+0.3*tr)
    s12=Ds*(1-0.3)*e12
    return s11,s22,s12,e11,e22
print("Far behind wake (x1=-18) profiles across x2:")
i=np.searchsorted(xs,-18.0)
for j in range(0,nyn+1,10):
    j=min(j,nyn-1)
    s11,s22,s12,e11,e22=stress(i,j)
    print(f" x2={ys[j]: .3f} s11={s11: .4f} s22={s22: .4f} s12={s12: .4f} u1={U1[i,j]: .4f} u2={U2[i,j]: .4f}")
print()
print("Just behind tip x1=-0.1:")
i=np.searchsorted(xs,-0.1)
for j in range(0,nyn+1,20):
    j=min(j,nyn-1)
    s11,s22,s12,e11,e22=stress(i,j)
    print(f" x2={ys[j]: .3f} s11={s11: .4f} s22={s22: .4f} s12={s12: .4f} u1={U1[i,j]: .4f} u2={U2[i,j]: .4f}")