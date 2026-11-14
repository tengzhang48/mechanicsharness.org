import numpy as np
from scipy.optimize import minimize_scalar
import json

# Recompute the solution and verify it independently
nu=0.3; q=0.5
K1,K2,K3=1.0,0.6,0.8
w=1-nu**2; v=1+nu
G0=w*(K1**2+K2**2)+v*K3**2
a1,a2,a3=0.5,-0.5,-0.5
B23,B32=-1.0,1.0
Cadv=2*(w*(a1*K1**2+a2*K2**2)+v*a3*K3**2)/G0
Ctilt=2*K2*K3*(w*B23+v*B32)/G0
Delta=Cadv**2+Ctilt**2

def fcf(s):
    th=2*np.pi*s
    D=1-2*q*np.cos(th)+q*q
    A=np.arctan2(q*np.sin(th),1-q*np.cos(th))
    return (0.5*Cadv*np.log(D)+Ctilt*A)/(np.pi*Delta)
def fpcf(s):
    th=2*np.pi*s
    D=1-2*q*np.cos(th)+q*q
    return (2*q/(Delta*D))*(-Cadv*np.sin(th)+Ctilt*(np.cos(th)-q))
def gcf(s):
    th=2*np.pi*s
    return (1-q*q)/(1-2*q*np.cos(th)+q*q)-1
def Hfcf(s):
    return (Ctilt*fpcf(s)-gcf(s))/Cadv

# extrema of f on one period
resmax=minimize_scalar(lambda s:-fcf(s), bounds=(-0.5,0.5), method='bounded', options={'xatol':1e-14})
resmin=minimize_scalar(lambda s: fcf(s), bounds=(-0.5,0.5), method='bounded', options={'xatol':1e-14})
print('f max',-resmax.fun,'at',resmax.x)
print('f min',resmin.fun,'at',resmin.x)
print('peak-to-peak',-resmax.fun-resmin.fun)

# verify kappa definitions from the Fourier series and the closed forms
N=8000
n=np.arange(1,N+1,dtype=float); qp=q**n
fc=(-Cadv)*qp/(np.pi*n*Delta); fs=(Ctilt)*qp/(np.pi*n*Delta)
def kappas_sum(s):
    th=2*np.pi*s
    H=np.sum(2*np.pi*n*(fc*np.cos(n*th)+fs*np.sin(n*th)))
    Fp=np.sum(2*np.pi*n*(-fc*np.sin(n*th)+fs*np.cos(n*th)))
    return np.array([-a1*(K1/K1)*H,
                     -a2*(K2/K1)*H+B23*(K3/K1)*Fp,
                     -a3*(K3/K1)*H+B32*(K2/K1)*Fp])
def kappas_cf(s):
    H,Hp=Hfcf(s),fpcf(s)
    return np.array([-a1*(K1/K1)*H,
                     -a2*(K2/K1)*H+B23*(K3/K1)*Hp,
                     -a3*(K3/K1)*H+B32*(K2/K1)*Hp])
s16=np.arange(16)/16
print('max |kappa sum - closed|', max(np.max(np.abs(kappas_sum(s)-kappas_cf(s))) for s in s16))

# load results and verify G balance
with open('results.json') as f: R=json.load(f)
print('json keys',R.keys())
print('PartA',R['partA'])
print('front_peak_to_peak',R['partB']['front_peak_to_peak'])
print('kappa3_max',R['partB']['kappa3_max'],'s',R['partB']['s_kappa3_max'])
for i,s in enumerate(s16):
    k=R['partB']['kappa1'][i]*K1+R['partB']['kappa2'][i]*K2+R['partB']['kappa3'][i]*K3
    dG=(2*w*(K1*R['partB']['kappa1'][i]+K2*R['partB']['kappa2'][i])+2*v*K3*R['partB']['kappa3'][i])/G0
    if abs(dG-gcf(s))>1e-10:
        print('G mismatch',s,dG,gcf(s))
print('G balance verified at all 16 points')
