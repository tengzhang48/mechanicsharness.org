import math,numpy as np,scipy.sparse.linalg as sla
import fem_solve as b, fem_quadratic as q, extract_k as ex
for L,nr in [(20,100),(40,106)]:
 m=q.elevate(b.build_mesh(nr,32,1e-4,L,h=1,H=1))
 K=q.assemble(m,E1=1,E2=1)
 F=q.consistent_loads(m,L,h=1,H=1)[:,1]
 fixed=np.array([2*m.right_tip,2*m.right_tip+1,2*m.right_top_corner])
 free=np.setdiff1d(np.arange(2*len(m.xy)),fixed)
 S=1/np.sqrt(K[free][:,free].diagonal());Kf=K[free][:,free]
 Ks=Kf.multiply(S[:,None]).multiply(S[None,:])
 u=np.zeros_like(F);u[free]=S*sla.spsolve(Ks,F[free]*S)
 ut=u.reshape(-1,2)[m.crack_top];ub=u.reshape(-1,2)[m.crack_bottom]
 delta=(ut-ub)[:,1]+1j*(ut-ub)[:,0]
 barE=1/(1-1/9);D=16/barE
 r=m.crack_r;ke=delta*2/D*np.sqrt(2*np.pi/r)
 k0,er,n=ex.extrapolate(r,ke,0,.02,2,.001)
 print(L,'amp0',abs(k0),'target',math.sqrt(12),'phase',np.angle(k0,deg=True))
 for rr in [.001,.005,.01,.02]:
  i=np.argmin(abs(r-rr));print(r[i],abs(ke[i]))
