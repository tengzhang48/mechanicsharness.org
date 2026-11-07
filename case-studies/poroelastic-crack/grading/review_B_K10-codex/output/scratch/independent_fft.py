"""Finite-crack resolvent check built from the independently derived symmetric symbol.
Unit L=G=c=1; unknowns: opening w and -p_y/sqrt(lambda).
This code imports no submitted solver. Saves boundary arrays for inspection.
"""
import argparse
import json
import os
import time
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg


def solve_grid(n,box,qs,tol=5e-11):
    nu=.2; nuu=.4; eta=(nuu-nu)/(1-nuu); C=1/(2*(1-nu))
    dx=box/n
    x=(np.arange(n)-n/2+.5)*dx
    ids=np.flatnonzero(abs(x)<1)
    nc=len(ids)
    k=2*np.pi*np.fft.rfftfreq(n,dx)

    def restriction(a,b,d):
        def mv(v):
            fields=np.zeros((2,n));fields[:,ids]=v.reshape(2,nc)
            hats=np.fft.rfft(fields,axis=1)
            out=np.vstack((a*hats[0]+b*hats[1],b*hats[0]+d*hats[1]))
            return np.fft.irfft(out,n=n,axis=1)[:,ids].ravel()
        return LinearOperator((2*nc,2*nc),matvec=mv,dtype=float)

    rhs=np.r_[np.ones(nc),np.zeros(nc)]
    aref=C*(1+eta)*k
    ref,info=cg(restriction(aref,0*k,np.ones_like(k)),rhs,
                M=restriction(1/np.maximum(aref,C/box),0*k,np.ones_like(k)),
                rtol=tol,atol=0,maxiter=3000)
    assert info==0,info
    wref=ref[:nc]
    padded_ref=np.zeros(n);padded_ref[ids]=wref
    stress_ref=-np.fft.irfft(aref*np.fft.rfft(padded_ref),n=n)
    outside=ids[-1]+1
    vref=dx*wref.sum()
    print(json.dumps(dict(kind='grid',N=n,box=box,dx=dx,crack_points=nc,
                          V0_grid=vref,V0_exact=np.pi*(1-nuu))),flush=True)
    previous=ref
    rows=[]
    for q in qs:
        started=time.monotonic()
        m=np.hypot(k,q)
        D=1+eta-2*eta*k*k/(m*(m+k))
        a=C*k*D
        b=np.sqrt(2*eta*C)*q*k/(m*(m+k))
        d=1/m
        M=restriction(a,b,d)
        det=a*d-b*b
        det[0]=C*d[0]/box
        ia=d/det;ib=-b/det;id_=a/det
        id_[0]=1/d[0]
        pre=restriction(ia,ib,id_)
        iterations=[0]
        def count(_):iterations[0]+=1
        sol,info=cg(M,rhs,x0=previous,M=pre,rtol=tol,atol=0,maxiter=3000,callback=count)
        relres=np.linalg.norm(M@sol-rhs)/np.linalg.norm(rhs)
        assert info==0,(info,relres)
        w=sol[:nc]; f=sol[nc:]
        V=dx*w.sum()
        # Elastic calibration has the same periodic box and endpoint quadrature.
        ratio=(1-nuu)/(1-nu)*w/wref
        tips={str(j):float(ratio[-j]) for j in [1,2,4,8,16,32]}
        padded=np.zeros((2,n));padded[:,ids]=sol.reshape(2,nc)
        hats=np.fft.rfft(padded,axis=1)
        stress=-np.fft.irfft(a*hats[0]+b*hats[1],n=n)
        stress_tips={str(j):float((stress[outside+j-1]+1)/(stress_ref[outside+j-1]+1)) for j in [1,2,4,8,16,32]}
        r=1-x[ids]
        # Fit the calibrated intensity at small r; finite-dx bias still remains.
        fitmask=(r/dx<17)
        Kfit=float(np.polynomial.polynomial.polyfit(r[fitmask],ratio[fitmask],2)[0])
        row=dict(kind='result',N=n,box=box,q=q,q_dx=q*dx,V=V,
                 A_same_grid=2*q*(V-vref)/np.sqrt(np.pi),
                 A_raw=2*q*(V-np.pi*(1-nuu))/np.sqrt(np.pi),
                 K_samples=tips,K_stress_samples=stress_tips,K_r0_fit=Kfit,relative_residual=relres,
                 iterations=iterations[0],seconds=time.monotonic()-started)
        rows.append(row);print(json.dumps(row),flush=True)
        previous=sol
        np.savez_compressed(f'scratch/fft_N{n}_box{box:g}_q{q:g}.npz',x=x[ids],w=w,flux_scaled=f,wref=wref,q=q,box=box,dx=dx)
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sizes',nargs='+',type=int,default=[32768,65536])
    p.add_argument('--box',type=float,default=32)
    p.add_argument('--qs',nargs='+',type=float,default=[20,80,320])
    args=p.parse_args()
    print(json.dumps(dict(pid=os.getpid(),started_utc=time.strftime('%H:%M:%S',time.gmtime()))),flush=True)
    rows=[solve_grid(n,args.box,args.qs) for n in args.sizes]
    if len(rows)>=2:
        ratio=args.sizes[-1]/args.sizes[-2]
        for j,q in enumerate(args.qs):
            A=(ratio*rows[-1][j]['A_same_grid']-rows[-2][j]['A_same_grid'])/(ratio-1)
            K=(ratio*rows[-1][j]['K_samples']['1']-rows[-2][j]['K_samples']['1'])/(ratio-1)
            print(json.dumps(dict(kind='dx_extrapolation',q=q,A=A,K_nearest=K)),flush=True)
