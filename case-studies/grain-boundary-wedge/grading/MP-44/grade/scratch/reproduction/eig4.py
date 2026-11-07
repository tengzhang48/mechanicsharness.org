import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

# Gauss for regular part
xg, wg = leggauss(60)

def Rreg(Z, X):
    return (X**2 - 4*Z*X - Z**2)/(Z+X)**3

def build_RM(N):
    edges = np.linspace(0,1,N+1)
    a = edges[:-1]; b = edges[1:]
    # collocation points: element midpoints
    Zc = 0.5*(a+b)
    # precompute per element regular integral as function of Z via Gauss nodes mapped
    Rmat = np.zeros((N,N)); Mmat=np.zeros((N,N))
    for i,Z in enumerate(Zc):
        # PV part per element: ln|(Z-a)/(Z-b)|
        pv = np.log(np.abs(Z-a)/np.abs(Z-b))
        # regular part per element
        for j in range(N):
            Xj = 0.5*(b[j]-a[j])*xg + 0.5*(a[j]+b[j])
            wr = 0.5*(b[j]-a[j])*wg
            reg = np.sum(wr*Rreg(Z,Xj))
            Rmat[i,j] = -(pv[j] + reg)
        # Build f from chi_j (indicator on element j) with BCs f(0)=0, f'(1)=f''(1)=0
        for j in range(N):
            # F0(Z)=int_0^Z chi_j = overlap of [0,Z] with [a_j,b_j]
            # We need F0,F1,F2 at Z. Compute by splitting.
            # element j contributes to F0(Z) if Z>a_j: len = clip(Z,a_j,b_j)-a_j
            def F0(Z, aj=a[j], bj=b[j]):
                return max(0.0, min(Z,bj)-aj)
            def F1(Z, aj=a[j], bj=b[j]):
                # int_0^Z F0(s) ds ; F0(s)=clip(s,a,b)-a for s in [a,b], 0 before a, (b-a) after b
                if Z<=aj: return 0.0
                zz=min(Z,bj)
                # int_a^zz (s-a) ds
                v=0.5*(zz-aj)**2
                if Z>bj: v += (b[j]-a[j])*(Z-bj)
                return v
            def F2(Z, aj=a[j], bj=b[j]):
                # int_0^Z F1 ds
                if Z<=aj: return 0.0
                zz=min(Z,bj)
                # int_a^zz 0.5(s-a)^2 ds = (1/6)(zz-a)^3
                v=(1.0/6.0)*(zz-aj)**3
                if Z>bj:
                    # F1 at bj = 0.5(b-a)^2 ; integrate from bj to Z: 0.5(b-a)^2 (Z-bj) + integral of extra term
                    L=bj-aj
                    v += 0.5*L**2*(Z-bj) + 0.5*L*(Z-bj)**2
                return v
            F0_1=F0(1.0); F1_1=F1(1.0)
            B = F0_1 - F1_1
            C = -F0_1/2.0
            Mmat[i,j] = F2(Z) + B*Z + C*Z**2
    return Rmat, Mmat, Zc

for N in [20,40,80,160,320]:
    Rmat,Mmat,Zc=build_RM(N)
    lam,vec=eig(Rmat,Mmat)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print("N=%3d lam1..4:"%N, np.round([c[0] for c in cand[:4]],6))