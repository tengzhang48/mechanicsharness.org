import numpy as np
import fem_j as FJ

ell=1.0
u,xs,ys,nxn,nyn,nid,shape,gauss=FJ.solve(ell,400,80,20.0,10.0)
U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
i0=np.searchsorted(xs,0.0)
print("Opening (top-bottom u2) behind tip:")
for i in range(i0-1,i0-10,-1):
    op=U2[i,nyn-1]-U2[i,0]
    print(f" x1={xs[i]: .4f}  open={op: .6e}")
print()
print("u2 at crack face x2=0 (x1<0), values:")
for i in range(i0-1,i0-8,-1):
    print(f" x1={xs[i]: .4f}  u2_bot={U2[i,0]: .6e} u2_top={U2[i,nyn-1]: .6e}")
print()
print("Far ahead stress check: compute sigma22 on x2=0 for x1>0")
# gradient at bottom row
dy=ys[1]-ys[0]
for i in range(i0,i0+6):
    u2y=(U2[i,1]-U2[i,0])/dy
    u1x=(U1[min(i+1,nxn-1),0]-U1[max(i-1,0),0])/(xs[min(i+1,nxn-1)]-xs[max(i-1,0)])
    e11=u1x; e22=u2y; tr=e11+e22
    Ds=1.0/(1-0.3**2)
    s22=1.0+Ds*((1-0.3)*e22+0.3*tr)
    print(f" x1={xs[i]: .4f} sigma22={s22: .6f}")