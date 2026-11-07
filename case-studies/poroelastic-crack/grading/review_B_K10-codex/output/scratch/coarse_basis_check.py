"""Audit the order of basis-resolution and small-time limits in MP-38."""
import os,sys,json,time
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,'scratch/MP-38')
from solve import spectral
print(json.dumps({'pid':os.getpid(),'started_utc':time.strftime('%H:%M:%S',time.gmtime())}),flush=True)
for n in [4,8]:
    rows=spectral(n,100000,.25,[80,320,1280,5120,20480])
    print(json.dumps({'basis_count':n,'rows_q_V_A_K':rows.tolist()}),flush=True)
