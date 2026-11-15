import numpy as np, fem, post, time
t0=time.time()
m=fem.build_mesh(x0=1e-6,q=1.15,Lx=3.0)
K=fem.assemble(m); F=fem.apply_loads(m,1.0,0.0); U=fem.solve(m,K,F)
r1,r2=post.interface_data(m,U)
np.savez('fineA_full.npz',coords=m['coords'],elems=m['elems'],mat=m['mat'],U=U,r1=r1,r2=r2)
print("done",time.time()-t0,flush=True)
