import numpy as np

nu=0.3; q=0.5
K1,K2,K3=1.0,0.6,0.8
w=1-nu**2; v=1+nu
G0=w*(K1**2+K2**2)+v*K3**2
a1,a2,a3=0.5,-0.5,-0.5
B23,B32=-1.0,1.0
Cadv=2*(w*(a1*K1**2+a2*K2**2)+v*a3*K3**2)/G0
Ctilt=2*K2*K3*(w*B23+v*B32)/G0
delta=Cadv**2+Ctilt**2
print("Cadv Ctilt",Cadv,Ctilt)

N=500
n=np.arange(1,N+1)
qpow=q**n
fc=(-Cadv)*qpow/(np.pi*n*delta)
fs=(Ctilt)*qpow/(np.pi*n*delta)
# Hf coefficients: multiply by 2*pi*n
hc=2*np.pi*n*fc
hs=2*np.pi*n*fs

def fsum(s):
    th=2*np.pi*s
    return np.sum(fc*np.cos(n*th)+fs*np.sin(n*th))
def fps(s):
    th=2*np.pi*s
    return np.sum(-2*np.pi*n*fc*np.sin(n*th)+2*np.pi*n*fs*np.cos(n*th))
def Hsum(s):
    th=2*np.pi*s
    return np.sum(hc*np.cos(n*th)+hs*np.sin(n*th))
def gsum(s):
    th=2*np.pi*s
    return np.sum(2*qpow*np.cos(n*th))

# closed forms
def fcf(s):
    th=2*np.pi*s
    D=1-2*q*np.cos(th)+q*q
    return (Cadv*0.5*np.log(D)+Ctilt*np.arctan2(q*np.sin(th),1-q*np.cos(th)))/(np.pi*delta)
def fpcf(s):
    th=2*np.pi*s
    D=1-2*q*np.cos(th)+q*q
    return (2*q/(np.pi*delta*D))*(-Cadv*np.sin(th)+Ctilt*(np.cos(th)-q))
def Hfcf(s):
    return (Ctilt*fpcf(s)-gsum(s))/Cadv

s=np.arange(16)/16
print("s   fsum    fcf    fps    fpcf   Hsum   Hfcf   g")
for ss in s:
    print(f"{ss:.6f} {fsum(ss): .8f} {fcf(ss): .8f} {fps(ss): .8f} {fpcf(ss): .8f} {Hsum(ss): .8f} {Hfcf(ss): .8f} {gsum(ss): .8f}")
print("max diff f", np.max(np.abs(fsum(s)-fcf(s))))
print("max diff Hf", np.max(np.abs(Hsum(s)-Hfcf(s))))
