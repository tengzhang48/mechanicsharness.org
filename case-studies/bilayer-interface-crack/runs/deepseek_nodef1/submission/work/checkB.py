import numpy as np, jdomain, time
t0=time.time()
d=np.load('betterB.npz')
coords,elems,mat,U=d['coords'],d['elems'],d['mat'],d['U']
for r1,r2 in [(0.02,0.1),(0.05,0.2)]:
    print(r1,r2,jdomain.jdomain(coords,elems,mat,U,r1,r2),flush=True)
print("time",time.time()-t0,flush=True)
