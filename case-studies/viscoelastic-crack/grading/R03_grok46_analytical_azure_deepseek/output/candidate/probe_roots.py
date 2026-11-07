import numpy as np, solve

def res2(V,a,eps,np_,ny):
    f,K,T0,W,S2,X,Y,K0 = solve.residuals(V,a,eps,npanel=np_,ny=ny)
    return f-0.1, K

def newton_safe(V,a,eps,np_,ny,itmax=14,verbose=True):
    for it in range(itmax):
        r0=np.array(res2(V,a,eps,np_,ny))
        if np.linalg.norm(r0)<1e-9: break
        da=1e-3*max(1,abs(a)); de=1e-3*max(1,abs(eps))
        r1=np.array(res2(V,a+da,eps,np_,ny))
        r2=np.array(res2(V,a,eps+de,np_,ny))
        J=np.column_stack([(r1-r0)/da,(r2-r0)/de])
        try: step=np.linalg.solve(J,-r0)
        except np.linalg.LinAlgError: break
        # damping: cap step
        nrm=np.linalg.norm(step)
        if nrm>0.5: step*=0.5/nrm
        a+=step[0]; eps+=step[1]
        if verbose: print('   it%d a=%.7f eps=%.7f |r|=%.2e'%(it,a,eps,np.linalg.norm(r0)))
    return a,eps,np.array(res2(V,a,eps,np_,ny))

if __name__=="__main__":
    for V in [1e-6, 0.001]:
        for np_,ny in [((140,90,140),80),((240,150,240),150),((360,230,360),230)]:
            a,eps,rr = newton_safe(V,0.441,0.521,np_,ny,verbose=False)
            print('V=%-8g res=%s -> a=%.7f eps=%.7f res=%.2e,%.2e'%(V,np_,a,eps,rr[0],rr[1]),flush=True)
