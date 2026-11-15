import numpy as np, jdomain
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])
eps=-0.06453177620670408
def fitT(r,T):
    n=len(r); M=np.zeros((2*n,2)); rhs=np.zeros(2*n)
    v=r**(-0.5+1j*eps); M[0:n,0]=np.real(v); M[0:n,1]=-np.imag(v); M[n:2*n,0]=np.imag(v); M[n:2*n,1]=np.real(v)
    rhs[0:n]=np.real(T); rhs[n:2*n]=np.imag(T)
    sol=np.linalg.lstsq(M,rhs,rcond=None)[0]
    return sol[0]+1j*sol[1],np.linalg.norm(M@sol-rhs)/np.linalg.norm(rhs)
for name in ['A','B']:
    d=np.load('better%s.npz'%name)
    print("="*20,"case",name)
    print("J:",[jdomain.jdomain(d['coords'],d['elems'],d['mat'],d['U'],r1,r2)[0] for r1,r2 in [(0.02,0.1),(0.05,0.2),(0.1,0.3)]])
    r1,r2=d['r1'],d['r2']; n=min(len(r1),len(r2)); r=r1[:n,0]
    T1=r1[:n,2]+1j*r1[:n,3]; T2=r2[:n,2]+1j*r2[:n,3]
    for tn,T in [('T1',T1),('T2',T2),('Tavg',0.5*(T1+T2))]:
        print(tn)
        for lo,hi in [(0.0005,0.005),(0.001,0.01),(0.002,0.02),(0.005,0.05),(0.01,0.1)]:
            sel=(r>=lo)&(r<=hi); D,res=fitT(r[sel],T[sel])
            print("  %.4g-%.4g arg %.2f abs %.3f res %.3f ph_h %.2f"%(lo,hi,np.degrees(np.angle(D)),abs(D),res,np.degrees(np.angle(D))+eps*np.log(0.5)*180/np.pi))
    # jump ratio to fine sub
    s=np.load('sub_jump_fine.npz'); rs,dxs,dys=s['r'],s['dx'],s['dy']; o=np.argsort(rs); rs=rs[o];dxs=dxs[o];dys=dys[o]
    ra,dxa,dya=d['r'],d['dx'],d['dy']; o=np.argsort(ra); ra=ra[o];dxa=dxa[o];dya=dya[o]
    sel=(rs>=3e-5)&(rs<=1e-3); rr=rs[sel]
    DX=np.interp(rr,ra,dxa); DY=np.interp(rr,ra,dya); DXs=np.interp(rr,rs,dxs); DYs=np.interp(rr,rs,dys)
    ratio=(DY+1j*DX)/(DYs+1j*DXs)
    print("jump ratio mean phase %.2f std %.2f mag %.3f  -> argKsub+ratio %.2f ph_h %.2f"%(np.mean(np.degrees(np.angle(ratio))),np.std(np.degrees(np.angle(ratio))),np.mean(np.abs(ratio)),np.degrees(np.angle(KA))+np.mean(np.degrees(np.angle(ratio))),np.degrees(np.angle(KA))+np.mean(np.degrees(np.angle(ratio)))+eps*np.log(0.5)*180/np.pi))
