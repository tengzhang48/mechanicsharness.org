"""Direct finite-crack boundary solve in the Laplace domain; no factorization.

python3 finite_crack_check.py (about several minutes; prints all report tables).
Optional --sizes, --box and --qs permit further refinement.  No files written.
The FFT applies the exact decaying-half-space Airy/diffusion surface operator.
Opening and normal pressure gradient are supported only on the slit.  See
REPORT.md for the derivation, calibration and extrapolation limitations.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import json
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from cauchy_check import face_calibration


def solve(n,box,qs):
    dx=box/n
    x=(np.arange(n)-n/2+.5)*dx
    slit=np.flatnonzero(np.abs(x)<1)
    nc=len(slit)
    k=2*np.pi*np.fft.rfftfreq(n,dx)
    C=1/(2*.8); eta=1/3

    def restricted(a,b,d):
        def action(v):
            full=np.zeros((2,n)); full[:,slit]=v.reshape(2,nc)
            f,g=np.fft.rfft(full,axis=1)
            out=np.fft.irfft(np.array([a*f+b*g,b*f+d*g]),n=n,axis=1)
            return out[:,slit].ravel()
        return LinearOperator((2*nc,2*nc),matvec=action,dtype=float)

    rhs=np.r_[np.ones(nc),np.zeros(nc)]
    refop=restricted(C*(1+eta)*k,0*k,np.ones_like(k))
    reference,info=cg(refop,rhs,M=restricted(1/np.maximum(C*(1+eta)*k,C/box),0*k,np.ones_like(k)),
                      rtol=2e-11,atol=0,maxiter=3000)
    assert info==0,info
    wref=reference[:nc]; vref=dx*wref.sum()
    print(json.dumps(dict(kind='elastic_calibration',N=n,box=box,dx=dx,
                          V0_grid=vref,V0_exact=.6*np.pi)),flush=True)
    previous=reference
    rows=[]
    for q in qs:
        m=np.hypot(k,q)
        a=C*k*(1+eta-2*eta*k*k/(m*(m+k)))
        b=np.sqrt(2*eta*C)*q*k/(m*(m+k))
        d=1/m
        op=restricted(a,b,d)
        det=a*d-b*b
        det[0]=C*d[0]/box
        ia=d/det; ib=-b/det; id_=a/det
        id_[0]=1/d[0]
        count=[0]
        def tick(_): count[0]+=1
        sol,info=cg(op,rhs,x0=previous,M=restricted(ia,ib,id_),
                    rtol=2e-11,atol=0,maxiter=4000,callback=tick)
        assert info==0,(n,q,info)
        w=sol[:nc]
        # Resolvent volume U(s)=s*Vhat(s): A=2*q*(U-U0)/sqrt(pi).
        A=2*q*dx*(w-wref).sum()/np.sqrt(np.pi)
        ratios=.75*w/wref
        # Calibrated opening ratios at the first cells approach the tip value.
        rr=1-x[slit][-12:]
        kfit=np.polynomial.polynomial.polyfit(rr,ratios[-12:],2)[0]
        row=dict(kind='direct',N=n,box=box,q=q,qdx=q*dx,A=float(A),
                 K_first=float(ratios[-1]),K_fit=float(kfit),iterations=count[0],
                 residual=float(np.linalg.norm(op@sol-rhs)/np.linalg.norm(rhs)))
        print(json.dumps(row),flush=True)
        rows.append(row);previous=sol
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sizes',nargs='+',type=int,default=[524288,1048576,2097152])
    p.add_argument('--box',type=float,default=32)
    p.add_argument('--qs',nargs='+',type=float,default=[80,160,320,640,1280])
    args=p.parse_args()
    face=face_calibration()  # check trusted flat-face term before the tip solve
    table=[solve(n,args.box,args.qs) for n in args.sizes]
    if len(table)>=2:
        ratio=args.sizes[-1]/args.sizes[-2]
        aa=[];kk=[]
        for j,q in enumerate(args.qs):
            a=(ratio*table[-1][j]['A']-table[-2][j]['A'])/(ratio-1)
            kval=(ratio*table[-1][j]['K_first']-table[-2][j]['K_first'])/(ratio-1)
            aa.append(a); kk.append(kval)
            print(json.dumps(dict(kind='dx_extrapolation',q=q,A=a,K=kval,tip=a-face)),flush=True)
        # The measured leading approach is q^(-1/2), not q^(-1).
        # This is an empirical extrapolation, with model sensitivity printed.
        q=np.array(args.qs,dtype=float)
        for degree in (1,2,3):
            for which in ('all','last4'):
                sl=slice(None) if which=='all' else slice(-4,None)
                design=np.array([q**(-j/2) for j in range(degree+1)]).T
                fit=np.linalg.lstsq(design[sl],np.array(aa)[sl],rcond=None)[0]
                fitk=np.linalg.lstsq(design[sl],np.array(kk)[sl],rcond=None)[0]
                print(json.dumps(dict(kind='large_q_fit',points=which,degree=degree,
                                      A=fit[0],tip=fit[0]-face,K=fitk[0],
                                      max_A_residual=float(np.max(np.abs(design[sl]@fit-np.array(aa)[sl]))))),flush=True)
