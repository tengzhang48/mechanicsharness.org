import numpy as np
import json

nu = 0.3
q = 0.5
K1 = 1.0
K2 = 0.6
K3 = 0.8
w = (1-nu**2)
v = (1+nu)
G0 = w*(K1**2 + K2**2) + v*K3**2

def sums(s, q=0.5, N=2000):
    """Return g, f, f', Hf using analytic formulas."""
    th = 2*np.pi*s
    D = 1 - 2*q*np.cos(th) + q*q
    g = (1-q*q)/D - 1
    # We'll compute f and Hf after coefficients known
    return th, D, g

def make_solution(a1,a2,a3,B, kind='IIIII'):
    Cadv = 2*(w*(a1*K1**2+a2*K2**2)+v*a3*K3**2)/G0
    if kind=='IIIII':
        # B[1,2]=B_II,III, B[2,1]=B_III,II in 1-indexed modes (2<->3)
        B23, B32 = B
        Ctilt = 2*K2*K3*(w*B23 + v*B32)/G0
    elif kind=='III':
        B12, B21 = B
        Ctilt = 2*K1*K2*(w*B12 + w*B21)/G0
    delta = Cadv**2 + Ctilt**2
    def f(s):
        th,D,g = sums(s)
        logterm = 0.5*np.log(D)
        atan = np.arctan2(q*np.sin(th), 1-q*np.cos(th))
        return (Cadv*logterm + Ctilt*atan)/(np.pi*delta)
    def fp(s):
        th,D,g = sums(s)
        return (2*q/(np.pi*delta*D))*(-Cadv*np.sin(th) + Ctilt*(np.cos(th)-q))
    def Hf(s):
        return (Ctilt*fp(s) - sums(s)[2])/Cadv
    def kappa(s):
        # kappa_alpha = (K_alpha/K1)[ -a_alpha Hf ] + sum_beta B_alpha,beta (K_beta/K1) fp
        kap = []
        Fp = fp(s); H = Hf(s)
        Ks=[K1,K2,K3]
        as_=[a1,a2,a3]
        for al in range(3):
            val = -(as_[al]*Ks[al]/K1)*H
            if kind=='IIIII':
                B23,B32=B
                if al==1: val += B23*(K3/K1)*Fp
                if al==2: val += B32*(K2/K1)*Fp
            elif kind=='III':
                B12,B21=B
                if al==0: val += B12*(K2/K1)*Fp
                if al==1: val += B21*(K1/K1)*Fp
            kap.append(val)
        return np.array(kap)
    return dict(Cadv=Cadv,Ctilt=Ctilt,delta=delta,f=f,fp=fp,Hf=Hf,kappa=kappa)

# Part A coefficients from same model
def partA(a1,a2,a3,B,kind='IIIII'):
    C1=-2*np.pi*a1*(K1/K1); C2=-2*np.pi*a2*(K2/K1); C3=-2*np.pi*a3*(K3/K1)
    if kind=='IIIII':
        B23,B32=B
        S1=0.0; S2=-2*np.pi*B23*(K3/K1); S3=-2*np.pi*B32*(K2/K1)
    elif kind=='III':
        B12,B21=B
        S1=-2*np.pi*B12*(K2/K1); S2=-2*np.pi*B21*(K1/K1); S3=0.0
    sol=make_solution(a1,a2,a3,B,kind)
    CG=-2*np.pi*sol['Cadv']; SG=-2*np.pi*sol['Ctilt']
    return dict(C1=C1,S1=S1,C2=C2,S2=S2,C3=C3,S3=S3,CG=CG,SG=SG),sol

cases = {
 'A_Rice_II_III_unstable': (0.5,-0.5,-0.5,(-1.0,1.0),'IIIII'),
 'B_I_II_stable_III_unstable': (0.5,0.5,-0.5,(-1.0,1.0),'IIIII'),
 'C_all_stable': (0.5,0.5,0.5,(-1.0,1.0),'IIIII'),
 'D_all_unstable': (-0.5,-0.5,-0.5,(-1.0,1.0),'IIIII'),
 'E_Rice_II_III_unstable_III_stable': (0.5,-0.5,0.5,(-1.0,1.0),'IIIII'),
 'F_CotterellRice_I_II_tilt': (0.5,0.5,-0.5,(1.5,-0.5),'III'),
}
s16 = np.arange(16)/16.0
for name,(a1,a2,a3,B,kind) in cases.items():
    pa,sol=partA(a1,a2,a3,B,kind)
    fs=sol['f'](s16)
    kap=sol['kappa'](s16)
    # sample fine for peak-to-peak and kappa3 max
    sf=np.linspace(-0.5,0.5,200001)
    ff=sol['f'](sf)
    k3=sol['kappa'](sf)[2]
    imax=np.argmax(k3)
    print('\n===',name)
    print('PartA', {k:round(float(v),6) for k,v in pa.items()})
    print('Cadv,Ctilt,delta',sol['Cadv'],sol['Ctilt'],sol['delta'])
    print('K10_scaled',np.sqrt((1-nu**2)/G0))
    print('front_peak_to_peak',ff.max()-ff.min(), ' fmin,fmax',ff.min(),ff.max(),'f(0)',sol['f'](0.0),'f(.5)',sol['f'](0.5))
    print('kappa3_max',k3[imax],'s=',sf[imax])
    print('f16',np.round(fs,6))
    print('k1',np.round(kap[0],6))
    print('k2',np.round(kap[1],6))
    print('k3',np.round(kap[2],6))
