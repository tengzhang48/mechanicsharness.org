import numpy as np, json
nu=0.3; q=0.5
K1,K2,K3=1.,0.6,0.8
w=1-nu**2; v=1+nu; G0=w*(K1**2+K2**2)+v*K3**2
a1,a2,a3=0.5,-0.5,-0.5; B23,B32=-1.,1.
Cadv=2*(w*(a1*K1**2+a2*K2**2)+v*a3*K3**2)/G0
Ctilt=2*K2*K3*(w*B23+v*B32)/G0
Delta=Cadv**2+Ctilt**2

def fcf(s):
    th=2*np.pi*s; D=1-2*q*np.cos(th)+q*q
    return (0.5*Cadv*np.log(D)+Ctilt*np.arctan2(q*np.sin(th),1-q*np.cos(th)))/(np.pi*Delta)
def fpcf(s):
    th=2*np.pi*s; D=1-2*q*np.cos(th)+q*q
    return (2*q/(Delta*D))*(Cadv*np.sin(th)+Ctilt*(np.cos(th)-q))
def gcf(s):
    th=2*np.pi*s
    return (1-q*q)/(1-2*q*np.cos(th)+q*q)-1
def Hfcf(s):
    return (Ctilt*fpcf(s)-gcf(s))/Cadv
def kapcf(s):
    H=Hfcf(s); Fp=fpcf(s)
    return np.array([-a1*(K1/K1)*H, -a2*(K2/K1)*H+B23*(K3/K1)*Fp, -a3*(K3/K1)*H+B32*(K2/K1)*Fp])

R=json.load(open('results.json'))
s16=np.arange(16)/16
maxdiff=0
for i,s in enumerate(s16):
    cf=kapcf(s)
    js=np.array([R['partB']['kappa1'][i],R['partB']['kappa2'][i],R['partB']['kappa3'][i]])
    maxdiff=max(maxdiff,np.max(np.abs(cf-js)))
    fd=abs(fcf(s)-R['partB']['f'][i])
    if fd>1e-8 or np.max(np.abs(cf-js))>1e-6:
        print('mismatch',s,'f',fcf(s),R['partB']['f'][i],'kap cf',cf,'json',js)
print('max kappa diff closed vs json',maxdiff)
print('f(0), f(1/16)',fcf(0),fcf(1/16))
# G balance
for i,s in enumerate(s16):
    dG=(2*w*(K1*R['partB']['kappa1'][i]+K2*R['partB']['kappa2'][i])+2*v*K3*R['partB']['kappa3'][i])/G0
    assert abs(dG-gcf(s))<1e-9,(s,dG,gcf(s))
print('G balance ok')
