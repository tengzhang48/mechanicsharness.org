import os,json,math,time
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from wilson_factor import wilson
print(json.dumps({'pid':os.getpid(),'started_utc':time.strftime('%H:%M:%S',time.gmtime())}),flush=True)
for eta in [0,1e-2,1e-3,1e-4]:
    rows=[wilson(eta,count,.01) for count in [16384,65536,262144]]
    n11=np.array([1/3,-2,8/3])@np.array([row['n11'] for row in rows])
    n21=np.array([1/3,-2,8/3])@np.array([row['n21'] for row in rows])
    K=n11/(1+eta)
    out=dict(eta=eta,K=K,n11=n11,n21=n21)
    if eta:
        out.update(K_scaled_deficit=(1-K)/eta,n11_first_derivative=(n11-1)/eta,
                   n21_scaled=n21/math.sqrt(eta/2),A_over_d=4/((1+eta)*math.sqrt(math.pi))*(4+math.pi*math.sqrt(2/eta)*n21))
    print(json.dumps(out),flush=True)
for iterations in [8,30,60]:
    row=wilson(1/3,65536,.01,iterations)
    print(json.dumps(dict(kind='iteration_count_variation',maxiter=iterations,**row)),flush=True)
