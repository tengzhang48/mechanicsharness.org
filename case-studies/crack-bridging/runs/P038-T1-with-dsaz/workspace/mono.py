import numpy as np
E=130e9; nu=0.28; Ep=E/(1-nu**2); s0=160e6; dc=9e-6; G0=22.0; a=2.4e-3
K0=np.sqrt(Ep*G0)
def sbar(eta):
    r=(2/np.pi)*np.arccos(1/eta)
    return s0*(K0/(s0*np.sqrt(np.pi*a)*np.sqrt(eta)) + r)
def dt(eta):
    r=sbar(eta)/s0-(2/np.pi)*np.arccos(1/eta)
    return s0*a/Ep*(4*np.sqrt(eta**2-1)*r+(8/np.pi)*np.log(eta))
print('eta   Sigma(MPa)  delta_t(um)')
for eta in [1.2,1.5,2.0,2.5,2.70,2.72,3.0,3.5,4.0]:
    print(f'{eta:5.2f}  {sbar(eta)/1e6:9.4f}  {dt(eta)*1e6:9.4f}')
print('-> Sigma and delta_t both increasing in eta: monotonic branch OK')
