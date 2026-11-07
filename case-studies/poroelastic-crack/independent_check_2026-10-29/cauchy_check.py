"""Numerical Cauchy/Riemann--Hilbert solve; no archived solver imports.

Run python3 cauchy_check.py.  This is an independent implementation, not an
independent physical model: see REPORT.md for the shared symbol and matching
identities.  GMRES solves the linear inverse-factor equation, not Wilson's
nonlinear spectral-factor iteration.  Only small JSON diagnostics are printed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
from scipy.sparse.linalg import LinearOperator, gmres


def face_calibration():
    integral, error = quad(lambda y: erfc(y/2), 0, np.inf, epsabs=2e-13)
    # The 1-D diffusion integral is 2/sqrt(pi).  Reciprocity and the
    # induced pressure multiply it by 8*(nu_u-nu)/(1+eta).
    face = 8*.2/(4/3)*integral
    print(json.dumps(dict(kind='face_calibration',erfc_integral=integral,
                         quadrature_error=error,face=face)),flush=True)
    assert abs(face-1.354055)<1e-6
    return face


def factor(count, scale, eta=1/3):
    theta=2*np.pi*np.arange(count)/count
    xi=scale*np.tan(theta/2)
    k=np.abs(xi); m=np.hypot(1,k)
    symbol=np.zeros((count,2,2),complex)
    symbol[:,0,0]=1+eta-2*eta*k*k/(m*(m+k))
    symbol[:,1,1]=1
    nz=k>0
    symbol[nz,0,1]=(np.sqrt(2*eta)*k[nz]/(m[nz]+k[nz])/
                    np.sqrt(1-1j*xi[nz])/np.sqrt(-1j*xi[nz]))
    symbol[:,1,0]=symbol[:,0,1].conj()
    symbol[count//2]=np.eye(2)
    delta=symbol-np.eye(2)
    modes=np.fft.fftfreq(count)<0
    modes[count//2]=False

    def cauchy_lower(v):
        coeff=np.fft.fft(v,axis=0)
        coeff[~modes]=0
        out=np.fft.ifft(coeff,axis=0)
        return out-out[count//2]  # zero at physical infinity

    def action(v):
        values=v.reshape(count,2)
        return (values+cauchy_lower(np.einsum('nij,nj->ni',delta,values))).ravel()

    operator=LinearOperator((2*count,2*count),matvec=action,dtype=complex)
    inverse=np.empty((count,2,2),complex)
    worst=0.
    for column in range(2):
        rhs=-cauchy_lower(delta[:,:,column]).ravel()
        solution,info=gmres(operator,rhs,rtol=2e-12,atol=1e-13,restart=25,maxiter=20)
        assert info==0, info
        worst=max(worst,float(np.linalg.norm(operator@solution-rhs)/np.linalg.norm(rhs)))
        inverse[:,:,column]=solution.reshape(count,2)+np.eye(2)[:,column]
    lower0=np.linalg.inv(inverse[0])
    n11=float(lower0[0,0].real); n21=float(lower0[1,0].real)
    row=dict(kind='cauchy',N=count,scale=scale,n11=n11,n21=n21,
             K=n11/(1+eta),tip=4*.2/((1+eta)*np.sqrt(np.pi))*np.pi*np.sqrt(2/eta)*n21,
             equation_residual=worst,energy_residual=n11*n11+n21*n21-(1+eta))
    print(json.dumps(row),flush=True)
    return row


if __name__=='__main__':
    face=face_calibration()
    estimates=[]
    for scale in (.01,.1):
        rows=[factor(n,scale) for n in (4096,16384,65536)]
        for weights,label in ((np.array([0,-1,2]),'linear'),
                              (np.array([1/3,-2,8/3]),'quadratic')):
            k=float(weights@np.array([r['K'] for r in rows]))
            tip=float(weights@np.array([r['tip'] for r in rows]))
            print(json.dumps(dict(kind='extrapolation',scale=scale,fit=label,
                                  K=k,tip=tip,total=face+tip)),flush=True)
            estimates.append((k,tip))
    assert all(abs(k-.807593)<3e-5 for k,_ in estimates)
    # A second identity, conditional on the trusted intensity, is diagnostic.
    n11=(4/3)*.807593
    n21=np.sqrt(4/3-n11*n11)
    tip=4*.2/((4/3)*np.sqrt(np.pi))*np.pi*np.sqrt(6)*n21
    print(json.dumps(dict(kind='conditional_energy_identity',input_K=.807593,
                          n21=n21,tip=float(tip))),flush=True)
