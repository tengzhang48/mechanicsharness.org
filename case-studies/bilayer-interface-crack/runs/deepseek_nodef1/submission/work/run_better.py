import numpy as np, fem, post, time, sys
t0=time.time()
for name,P,M in [('A',1,0),('B',0,1)]:
    m=fem.build_mesh(x0=1e-6,q=1.07,Lx=3.0)
    print(name,"mesh",m['ntot'],m['elems'].shape[0],time.time()-t0,flush=True)
    K=fem.assemble(m); print(name,"assemble",time.time()-t0,flush=True)
    F=fem.apply_loads(m,P,M); U=fem.solve(m,K,F); print(name,"solve",time.time()-t0,flush=True)
    r,dx,dy=post.crack_jump(m,U); o=np.argsort(r)
    r1,r2=post.interface_data(m,U)
    np.savez('better%s.npz'%name,coords=m['coords'],elems=m['elems'],mat=m['mat'],U=U,r=r[o],dx=dx[o],dy=dy[o],r1=r1,r2=r2)
    print(name,"done",time.time()-t0,flush=True)
print("ALL DONE",flush=True)
