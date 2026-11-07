"""Independent right spectral factor by Wilson iteration.
No imports from supplied solvers; solves N = F^* F, F analytic below.
Fourier convention: negative circle modes are lower half-plane functions.
"""
import argparse
import json
import os
import time
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np


def wilson(eta=1/3, count=16384, scale=.03, maxiter=30):
    theta=2*np.pi*np.arange(count)/count
    xi=scale*np.tan(theta/2)
    k=np.abs(xi)
    m=np.hypot(1.,k)
    # Construct directly from the symmetrized displacement/flux symbol,
    # including its diagonal analytic factors.
    N=np.zeros((count,2,2),complex)
    N[:,0,0]=1+eta-2*eta*k*k/(m*(m+k))
    N[:,1,1]=1
    ids=k>0
    N[ids,0,1]=np.sqrt(2*eta)*k[ids]/(m[ids]+k[ids])/(np.sqrt(1-1j*xi[ids])*np.sqrt(-1j*xi[ids]))
    N[:,1,0]=N[:,0,1].conj()
    N[count//2]=np.eye(2)
    freq=np.fft.fftfreq(count)
    lower=freq<0
    lower[count//2]=False

    def half(values):
        modes=np.fft.fft(values,axis=0)
        constant=modes[0].copy()/2
        modes[~lower]=0
        # Wilson projection has half of the constant Fourier mode.
        modes[0]=constant
        return np.fft.ifft(modes,axis=0)

    F=np.broadcast_to(np.eye(2),(count,2,2)).astype(complex).copy()
    for iteration in range(maxiter):
        inv=np.linalg.inv(F)
        Q=inv.conj().transpose(0,2,1) @ N @ inv
        error=float(np.max(np.abs(Q-np.eye(2))))
        if error<1e-11:
            break
        F=half(Q+np.eye(2)) @ F
    # At physical infinity, N=I so F(infinity) is unitary; remove that unitary.
    F=np.linalg.inv(F[count//2]) @ F
    actual=float(np.max(np.abs(F.conj().transpose(0,2,1) @ F-N)))
    coeffs=np.fft.fft(F,axis=0)/count
    leakage=float(np.max(np.abs(coeffs[freq>0])))
    zero=F[0]
    n11=float(zero[0,0].real)
    n21=float(zero[1,0].real)
    d=.8*eta/(1+eta) # keep drained nu=.2 as coupling is varied
    A=4*d/((1+eta)*np.sqrt(np.pi))*(4+np.pi*np.sqrt(2/eta)*n21) if eta else 0
    return dict(method='Wilson right spectral factor', count=count,scale=scale,eta=eta,
                iterations=iteration,n11=n11,n21=n21,K=n11/(1+eta),A=A,
                max_factor_residual=actual,max_positive_mode=leakage,
                norm_residual=n11*n11+n21*n21-(1+eta),zero_imag=float(abs(zero.imag).max()))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--counts',nargs='+',type=int,default=[4096,16384,65536])
    p.add_argument('--scales',nargs='+',type=float,default=[.01,.1,1.])
    p.add_argument('--eta',type=float,default=1/3)
    p.add_argument('--iterations',type=int,default=30)
    args=p.parse_args()
    print(json.dumps({'pid':os.getpid(),'started_utc':time.strftime('%H:%M:%S',time.gmtime())}),flush=True)
    for scale in args.scales:
        rows=[]
        for count in args.counts:
            row=wilson(args.eta,count,scale,args.iterations)
            rows.append(row)
            print(json.dumps(row),flush=True)
        if len(rows)==3 and args.counts[1]==4*args.counts[0] and args.counts[2]==4*args.counts[1]:
            w=np.array([1/3,-2,8/3])
            n11=w @ np.array([r['n11'] for r in rows])
            n21=w @ np.array([r['n21'] for r in rows])
            print(json.dumps({'scale':scale,'quadratic_extrap_n11':n11,'quadratic_extrap_n21':n21,
                              'K':n11/(1+args.eta),'A':4*(.8*args.eta/(1+args.eta))/((1+args.eta)*np.sqrt(np.pi))*(4+np.pi*np.sqrt(2/args.eta)*n21),
                              'norm_residual':n11*n11+n21*n21-(1+args.eta)}),flush=True)
