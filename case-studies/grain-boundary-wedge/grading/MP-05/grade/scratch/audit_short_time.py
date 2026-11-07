import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent / 'reproduction'))
from solver10 import Res10

for n in [48, 80, 120]:
    r = Res10(n, Mquad=600, Kleg=200)
    f = r._frows([0.3, 0.5, 0.7, 1.0])
    print('N', n)
    for p in [1e3, 1e4, 1e5, 1e6, 1e7]:
        stress_hat = f @ r.solve(p)
        print('p', p, 'p*(p*sigmahat-1)', (p*(p*stress_hat-1)).tolist())
print('boundary_layer_limit', (16/np.array([0.3,0.5,0.7,1.0])**3).tolist())
