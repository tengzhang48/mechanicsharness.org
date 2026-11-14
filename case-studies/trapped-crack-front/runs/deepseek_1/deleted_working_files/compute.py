import numpy as np
import json, sys

nu = 0.3
E = 1.0  # arbitrary; K10_scaled uses E*Gamma0
q = 0.5
K1 = 1.0
K2 = 0.6
K3 = 0.8

w = (1-nu**2)/E       # weight for modes I, II in G
v = (1+nu)/E         # weight for mode III in G
G0 = w*(K1**2 + K2**2) + v*K3**2   # times K1^0^2; actually G0 = ...*K1^0^2/E

print("w,v,G0", w, v, G0)

# Model A: diagonal advance coefficients
# delta K_alpha = - a_alpha K_alpha^0 |k| delta a + tilt
# choose a1=1/2, a2=a3=-1/2 (Rice: mode I stable, II/III unstable)
def partA(a1,a2,a3,B23=-1.0,B32=1.0):
    C1 = -2*np.pi*a1*(K1/K1)
    C2 = -2*np.pi*a2*(K2/K1)
    C3 = -2*np.pi*a3*(K3/K1)
    S1 = 2*np.pi*(B23*0 + B32*0)  # no I-II or I-III tilt
    S2 = 2*np.pi*(B23*(K3/K1))   # delta K2 = B23 K3 a'
    S3 = 2*np.pi*(B32*(K2/K1))   # delta K3 = B32 K2 a'
    # G response
    Cadv = 2*(w*(a1*K1**2+a2*K2**2)+v*a3*K3**2)/G0
    Ctilt = 2*K2*K3*(v*B32 + w*B23)/G0  # check: G tilt = 2K2K3 a'(v B32 + w B23)? derive below
    # Derivation: delta G = 2w K2 delta K2 + 2v K3 delta K3
    # delta K2 = B23 K3 a', delta K3 = B32 K2 a'
    # = 2w K2 B23 K3 a' + 2v K3 B32 K2 a' = 2K2K3 a'(w B23 + v B32)
    CG = -2*np.pi*Cadv
    SG = -2*np.pi*Ctilt
    return dict(C1=C1,S1=S1,C2=C2,S2=S2,C3=C3,S3=S3,CG=CG,SG=SG,
                Cadv=Cadv,Ctilt=Ctilt)

for name,(a1,a2,a3) in {
    "rice_II_III_unstable": (0.5,-0.5,-0.5),
    "all_stable": (0.5,0.5,0.5),
    "all_unstable": (-0.5,-0.5,-0.5),
    "I_stable_II_stable_III_unstable": (0.5,0.5,-0.5),
}.items():
    pa=partA(a1,a2,a3)
    print("\n",name)
    for k in ["C1","S1","C2","S2","C3","S3","CG","SG","Cadv","Ctilt"]:
        print(f"  {k}={pa[k]:.6f}")
    print("  K10_scaled", np.sqrt((1-nu**2)/(G0)))  # = sqrt((1-nu^2)/(G0)) where G0 in units K1^2/E
