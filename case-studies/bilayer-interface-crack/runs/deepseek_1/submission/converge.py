import numpy as np, time
from fem_bilayer import Mesh, Solver
from analyze_case import eps,Estr,Ceff,fit_q,analyze_cod,analyze_stress,solve_case

def run(P,M,params,label):
    p,nx,beta,ny1,ny2,L,gamma=params
    t0=time.time()
    m,s,u=solve_case(P,M,p=p,nx=nx,beta_x=beta,ny1=ny1,ny2=ny2,L=L,gamma=gamma)
    r,d,Q,outc=analyze_cod(s,u)
    xs,Qs,outs=analyze_stress(s,u)
    # pick w2 deg2 for COD and stress
    Kc=outc['w2'][2][2]; Ks=outs['w2'][2][2] if 2 in outs['w2'] else None
    print('%s params %s dofs %d time %.1f'%(label,params,m.ndof,time.time()-t0))
    print('   COD w2 deg2 |K| %.5f psi %.3f ; stress w2 |K| %.5f psi %.3f'%(abs(Kc),np.degrees(np.angle(Kc)),abs(Ks),np.degrees(np.angle(Ks))))
    return Kc,Ks

if __name__=='__main__':
    for tag,P,M in [('A',1.0,0.0),('B',0.0,1.0)]:
        for params,label in [((2,100,3.0,30,36,10.0,2.0),'coarse-p2'),((3,60,2.8,24,30,12.0,2.0),'medium-p3')]:
            run(P,M,params,tag+' '+label)
