"""Williams-series extrapolations of stored crack-face opening data."""
import sys
import numpy as np
import bilayer_fe as b

files=sys.argv[1:] or ['run_coarse.npz','run_medium.npz']
for fn in files:
 d=np.load(fn);r=d['r'];K=d['Kline'];print('\n',fn)
 for family in [1,2,3]:
  for rmax in [.02,.04,.08,.16,.32,.64]:
   for rmin in [0,1e-4,3e-4,1e-3]:
    mask=(r<=rmax)&(r>rmin)
    rr=r[mask,None]
    if family==1:
     A=np.concatenate([np.ones_like(rr),rr**(1-b.LAMBDA)],axis=1)
    elif family==2:
     A=np.concatenate([np.ones_like(rr),rr**(1-b.LAMBDA),rr],axis=1)
    else:
     A=np.concatenate([np.ones_like(rr),rr**(1-b.LAMBDA),rr,rr**(2-b.LAMBDA),rr**2],axis=1)
    c=np.linalg.lstsq(A,K[mask],rcond=None)[0][0]
    print(f'f={family} max={rmax} min={rmin} n={mask.sum()}', ' '.join(f'{abs(v):.7f}@{np.degrees(np.angle(v)):.5f}' for v in c))
