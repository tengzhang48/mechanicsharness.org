import numpy as np, fem, post, subproblem, time, sys
t0=time.time()
m,U=subproblem.run(x0=1e-6,q=1.15,a=0.5)
r,dx,dy=post.crack_jump(m,U); o=np.argsort(r); np.savez('sub_jump_fine.npz',r=r[o],dx=dx[o],dy=dy[o])
print("sub done",time.time()-t0,len(r),flush=True)
for name,P,M in [('A',1,0),('B',0,1)]:
    m=fem.build_mesh(x0=1e-6,q=1.15,Lx=3.0)
    K=fem.assemble(m); F=fem.apply_loads(m,P,M); U=fem.solve(m,K,F)
    r,dx,dy=post.crack_jump(m,U); o=np.argsort(r); np.savez('act%s_jump_fine.npz'%name,r=r[o],dx=dx[o],dy=dy[o])
    print("act",name,"done",time.time()-t0,m['ntot'],len(r),flush=True)
print("ALL DONE",flush=True)
