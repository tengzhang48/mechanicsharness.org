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

N=2000
n=np.arange(1,N+1)
qp=q**n
fc=(-Cadv)*qp/(np.pi*n*delta)
fs=(Ctilt)*qp/(np.pi*n*delta)

def fsum(s):
    th=2*np.pi*s
    return np.sum(fc*np.cos(n*th)+fs*np.sin(n*th))
def fps(s):
    th=2*np.pi*s
    return np.sum(-2*np.pi*n*fc*np.sin(n*th)+2*np.pi*n*fs*np.cos(n*th))
def Hsum(s):
    th=2*np.pi*s
    return np.sum(2*np.pi*n*fc*np.cos(n*th)+2*np.pi*n*fs*np.sin(n*th))
def gsum(s):
    th=2*np.pi*s
    return np.sum(2*qp*np.cos(n*th))

def kap(s):
    f=fsum(s); fp=fps(s); H=Hsum(s); g=gsum(s)
    k1 = -a1*1.0*H
    k2 = -a2*(K2/K1)*H + B23*(K3/K1)*fp
    k3 = -a3*(K3/K1)*H + B32*(K2/K1)*fp
    return f,fp,H,g,np.array([k1,k2,k3])

s16=np.arange(16)/16
F=[]; K1v=[];K2v=[];K3v=[]
for s in s16:
    f,fp,H,g,k=kap(s)
    F.append(f);K1v.append(k[0]);K2v.append(k[1]);K3v.append(k[2])
    # verify k3 = (K3/K1)*(-a3)*H + ... and equilibrium
print("Cadv,Ctilt,delta",Cadv,Ctilt,delta)
print("f",np.array2string(np.array(F),precision=10,separator=', '))
print("k1",np.array2string(np.array(K1v),precision=10,separator=', '))
print("k2",np.array2string(np.array(K2v),precision=10,separator=', '))
print("k3",np.array2string(np.array(K3v),precision=10,separator=', '))
# fine grid for peak-to-peak and kappa3 max
sf=np.linspace(-0.5,0.5,400001)
ff=np.array([fsum(x) for x in sf])
k3f=np.array([kap(x)[4][2] for x in sf])
imax=np.argmax(k3f)
print("peak-to-peak",ff.max()-ff.min(),"min,max",ff.min(),ff.max())
print("kappa3_max",k3f[imax],"s",sf[imax])
print("k3 = 3.316*g? ratio", np.mean(k3f/(gsum(sf)+1e-30)))
# Part A
C1=-2*np.pi*a1
C2=-2*np.pi*a2*(K2/K1)
C3=-2*np.pi*a3*(K3/K1)
S1=0.0
S2=-2*np.pi*B23*(K3/K1)
S3=-2*np.pi*B32*(K2/K1)
CG=-2*np.pi*Cadv
SG=-2*np.pi*Ctilt
print("PartA",C1,S1,C2,S2,C3,S3,CG,SG)
print("K10_scaled",np.sqrt((1-nu**2)/G0))
