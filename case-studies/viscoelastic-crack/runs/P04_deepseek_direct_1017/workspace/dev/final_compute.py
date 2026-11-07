import numpy as np, json, time
from global_bvp_g import solve_global
from dissipation_g import fields_T, gamma_D_ratio
from dissipation_exp import gamma_D_ratio_exp
from dissipation_psi import gamma_D_psi

def top_row(r,V,lam=0.1):
    X=r['X']; Y=r['Y']; U=r['U']; Ny=len(Y)-1; a=r['a']
    hm=Y[Ny-1]-Y[Ny-2]; hp=Y[Ny]-Y[Ny-1]
    # 3-point one-sided derivative at Y=pi (points Y_{N-2},Y_{N-1},pi)
    cA=(hm+2*hp)/(hp*(hm+hp)); cB=-(hm+hp)/(hm*hp); cC=hp/(hm*(hm+hp))
    E2b=np.pi*(cA*a+cB*U[-1,:]+cC*U[-2,:])
    E1b=np.zeros_like(E2b)  # w_x=0 at top
    # ODE for T along X at top
    def ode(E,er):
        t=np.zeros_like(E); tp=er; ep=er
        for k in range(len(E)-1,-1,-1):
            i=k+1; h=X[i+1]-X[i]
            R=E[k]-V*(ep-E[k])/h
            t[k]=(R+lam*V*tp/h)/(1+lam*V/h)
            tp=t[k];ep=E[k]
        return t
    T1b=ode(E1b,0.0); T2b=ode(E2b,a)
    return E1b,E2b,T1b,T2b

def sample_fields(r,Xs,ys,V):
    X=r['X']; Y=r['Y']; Xi=X[1:-1]; Yi=Y[1:-1]
    E1,E2,T1,T2=fields_T(r,V)
    E1b,E2b,T1b,T2b=top_row(r,V)
    out=[]
    for xq in Xs:
        for yq in ys:
            Yq=np.pi*yq
            if yq>=1.0-1e-12:
                e1=0.0; e2=np.interp(xq,Xi,E2b); t1=0.0; t2=np.interp(xq,Xi,T2b)
            else:
                # bilinear interpolation
                def bil(F):
                    # interpolate along X for each Y row then along Y
                    fi=np.array([np.interp(xq,Xi,F[j,:]) for j in range(F.shape[0])])
                    return np.interp(Yq,Yi,fi)
                e1=bil(E1); e2=bil(E2); t1=bil(T1); t2=bil(T2)
            out.append(dict(X=float(xq),y_over_h=float(yq),T1=float(t1),T2=float(t2),E1=float(e1),E2=float(e2)))
    return out

def compute(V,a,eps,grid,xs=(-4.0,-2.0,0.0,2.0,4.0),ys=(0.1,0.5,1.0)):
    t0=time.time()
    r=solve_global(a,eps,V,**grid)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps));ue=U0[it+1]
    gpsi=gamma_D_psi(r,V) if V else 0.0
    gold=gamma_D_ratio(r,V) if V else 0.0
    gexp=gamma_D_ratio_exp(r,V) if V else 0.0
    fields=sample_fields(r,list(xs),list(ys),V)
    return dict(V=V,a=a,eps=eps,ue=ue,R=a*a/0.2,gpsi=gpsi,gold=gold,gexp=gexp,fields=fields,t=time.time()-t0)

if __name__=='__main__':
    grid=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
    cases=[(0.001,0.451575,0.548932),(0.1,0.5205,0.755)]
    res={}
    for V,a,eps in cases:
        d=compute(V,a,eps,grid); res[V]=d
        print('V',V,'a',a,'eps',eps,'ue %.6f R %.5f gpsi %.5f gold %.5f gexp %.5f t %.1f'%(d['ue'],d['R'],d['gpsi'],d['gold'],d['gexp'],d['t']))
        print('  energy checks: R-1-gpsi %.5f  rel %.4f'%(d['R']-1-d['gpsi'],(d['R']-1-d['gpsi'])/max(1,d['R'])))
    json.dump(res,open('/workspace/dev/final_raw.json','w'),indent=1,default=float)
