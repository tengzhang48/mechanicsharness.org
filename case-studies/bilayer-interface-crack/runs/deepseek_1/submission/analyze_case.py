import numpy as np, time, json
from fem_bilayer import Mesh, Solver

eps = (1/(2*np.pi))*np.log((5/3/(3.375)+1/0.375)/(5/3/(0.375)+1/3.375))
mu1=3.375; mu2=0.375; nu=1/3
E1s=mu1*2*(1+nu)/(1-nu**2); E2s=mu2*2*(1+nu)/(1-nu**2)
Estr=1/(0.5*(1/E1s+1/E2s))
beta=0.2
Ceff=8/(Estr*np.cosh(np.pi*eps))/(1+2j*eps)

def solve_case(P,M,p=3,nx=100,beta_x=3.1,ny1=36,ny2=44,L=10.0,gamma=2.0):
    m=Mesh(h=1.0,H=2.0,L=L,p=p,nxL=nx,nxR=nx,beta_x=beta_x,ny1=ny1,ny2=ny2,gamma_y=gamma,crack=True)
    s=Solver(m,E1=9.0,E2=1.0,nu1=nu,nu2=nu)
    s.assemble(); u=s.solve(P=P,M=M)
    return m,s,u

def fit_q(r,Q,rmin,rmax,deg):
    sel=(r>=rmin)&(r<=rmax); rr=r[sel]; qq=Q[sel]
    if len(rr)<deg+1: return np.nan
    A=np.vander(rr,deg+1,increasing=True)
    cr=np.linalg.lstsq(A,qq.real,rcond=None)[0][0]
    ci=np.linalg.lstsq(A,qq.imag,rcond=None)[0][0]
    return cr+1j*ci

def analyze_cod(s,u):
    r,d=s.cod(u)
    Q=d*np.sqrt(2*np.pi/r)/Ceff*r**(-1j*eps)
    out={}
    for tag,(a,b) in {'w1':(0.005,0.05),'w2':(0.01,0.1),'w3':(0.02,0.1),'w4':(0.005,0.1)}.items():
        out[tag]={}
        for deg in [1,2,3]:
            K0=fit_q(r,Q,a,b,deg)
            out[tag][deg]=(abs(K0),np.degrees(np.angle(K0)),K0)
    return r,d,Q,out

def analyze_stress(s,u):
    # choose x values log-spaced in a range, avoiding element boundaries perhaps
    xs=np.geomspace(0.002,0.2,30)
    vals=[]
    for x in xs:
        sig=s.stress_interface(u,x)
        vals.append(sig[1]+1j*sig[2])
    vals=np.array(vals)
    Qs=vals*np.sqrt(2*np.pi*xs)*xs**(-1j*eps)
    out={}
    for tag,(a,b) in {'w1':(0.005,0.05),'w2':(0.01,0.1),'w3':(0.02,0.1),'w4':(0.005,0.1)}.items():
        out[tag]={}
        for deg in [1,2]:
            K0=fit_q(xs,Qs,a,b,deg)
            out[tag][deg]=(abs(K0),np.degrees(np.angle(K0)),K0)
    return xs,Qs,out

if __name__=='__main__':
    import sys
    tag=sys.argv[1] if len(sys.argv)>1 else 'A'
    if tag=='A': P,M=1.0,0.0
    elif tag=='B': P,M=0.0,1.0
    elif tag=='C': P,M=1.0,1.0
    else: raise SystemExit('tag must be A, B or C')
    t0=time.time(); m,s,u=solve_case(P,M); print('solve time',time.time()-t0,'dofs',m.ndof)
    r,d,Q,outc=analyze_cod(s,u)
    print('COD fits (window: deg1 deg2 deg3 -> |K|,psi):')
    for w,v in outc.items():
        print(' ',w, ['%.5f/%.3f'%(v[d][0],v[d][1]) for d in [1,2,3]])
    xs,Qs,outs=analyze_stress(s,u)
    print('Stress fits:')
    for w,v in outs.items():
        print(' ',w, ['%.5f/%.3f'%(v[d][0],v[d][1]) for d in [1,2]])
    # print raw stress Q at selected x
    for x in [0.002,0.005,0.01,0.02,0.05,0.1,0.2]:
        k=np.argmin(abs(xs-x)); print(' x %.4f Qs %s |K| %.5f psi %.3f'%(xs[k],Qs[k],abs(Qs[k]),np.degrees(np.angle(Qs[k]))))
    np.savez('data/cod_%s.npz'%tag,r=r,d=d)
