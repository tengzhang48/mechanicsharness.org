import numpy as np
np.set_printoptions(precision=8, suppress=False, linewidth=200)
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])

lam = s+1
# basis functions of potentials and their derivatives
def pot_stress_terms(x,y):
    """Return arrays d(sxx)/dA... for 8 real params. Complex stress expressions linear in A,B,C,D and conjugates.
       We'll just return the complex coefficient of each of the 8 real params for F and trace."""
    z=x+1j*y; zc=np.conj(z)
    # We write everything as linear in A, Abar, B, Bbar, C, Cbar, D, Dbar.
    # Then convert to real params: A=a+ib etc.
    # Build coefficient vectors for F and trace for each of the 8 complex monomials.
    # F = Phi + Phic + z*Phipc + Psibar
    # Phi = lam*A*z^(lam-1) + lbar*B*z^(lbar-1)
    # Phic = conj(Phi) = lbar*Abar*zc^(lbar-1) + lam*Bbar*zc^(lam-1)
    # Phip = lam*(lam-1)*A*z^(lam-2) + lbar*(lbar-1)*B*z^(lbar-2)
    # Phipc = lbar*(lbar-1)*Abar*zc^(lbar-2) + lam*(lam-1)*Bbar*zc^(lam-2)
    # Psi = (lam+1)*lam*C*z^(lam-1) + (lbar+1)*lbar*D*z^(lbar-1)
    # Psibar = conj(Psi) = (lbar+1)*lbar*Cbar*zc^(lbar-1) + (lam+1)*lam*Dbar*zc^(lam-1)
    l=lam; lb=np.conj(lam)
    zl=z**(l-1); zlb=z**(lb-1); zl2=z**(l-2); zlb2=z**(lb-2)
    zcl=zc**(l-1); zclb=zc**(lb-1); zcl2=zc**(l-2); zclb2=zc**(lb-2)
    # coefficients for F in terms of A,Abar,B,Bbar,C,Cbar,D,Dbar
    cF = {
      'A': l*zl,
      'Abar': lb*zclb + z*lb*(lb-1)*zclb2,
      'B': lb*zlb,
      'Bbar': l*zcl + z*l*(l-1)*zcl2,
      'C': (l+1)*l*zl,
      'Cbar': 0,
      'D': (lb+1)*lb*zlb,
      'Dbar': 0,
    }
    # Wait: Psibar terms for Cbar,Dbar are nonzero. Fix:
    cF['Cbar'] = (lb+1)*lb*zclb
    cF['Dbar'] = (l+1)*l*zcl
    # trace = 2(Phi+Phic)
    cT = {
      'A': 2*l*zl,
      'Abar': 2*lb*zclb,
      'B': 2*lb*zlb,
      'Bbar': 2*l*zcl,
      'C': 0,'Cbar':0,'D':0,'Dbar':0,
    }
    return cF,cT

def build_system(cAiry, pts):
    rows=[]; rhs=[]
    names=['A','Abar','B','Bbar','C','Cbar','D','Dbar']
    for (x,y) in pts:
        r=np.hypot(x,y); th=np.arctan2(y,x)
        f=cAiry@basis(th,s); fp=cAiry@dbasis(th,s); fpp=cAiry@ddbasis(th,s)
        Prr=(s+2)*f+fpp; Ptt=(s+2)*(s+1)*f; Prt=-(s+1)*fp
        rs=r**s
        Sxx=rs*(Prr*np.cos(th)**2+Ptt*np.sin(th)**2-2*Prt*np.sin(th)*np.cos(th))
        Syy=rs*(Prr*np.sin(th)**2+Ptt*np.cos(th)**2+2*Prt*np.sin(th)*np.cos(th))
        Sxy=rs*((Prr-Ptt)*np.sin(th)*np.cos(th)+Prt*(np.cos(th)**2-np.sin(th)**2))
        target=[np.real(Sxx),np.real(Syy),np.real(Sxy)]
        cF,cT=pot_stress_terms(x,y)
        # F and T are complex linear combos.  Convert to real params:
        # for each complex name, write coeff * val; val is complex unknown.
        # We use 8 complex monomials but they are conjugates; treat as independent complex? No.
        # Instead use real params: A=Ar+iAi etc. For each term coeff*X, with X complex.
        # Build 3 x 8 real matrix directly by evaluating derivative wrt each real param.
        # We'll use the identity: coeff*X + coeff2*conj(X) for pairs.
        # Simpler: use complex monomial values as 8 independent real params? They are constrained.
        # Let's instead take A,B,C,D complex directly and solve complex-linear system? 
        # F = cA*A + cAb*Abar + ... is not complex-linear but is real-linear.
        # Use 8 real unknowns: Ar,Ai,Br,Bi,Cr,Ci,Dr,Di.
        # Map: A=Ar+iAi, Abar=Ar-iAi, etc.
        M=np.zeros((3,8))
        for j in range(8):
            # perturb real param j by 1
            vals={'A':0j,'Abar':0j,'B':0j,'Bbar':0j,'C':0j,'Cbar':0j,'D':0j,'Dbar':0j}
            keys=['A','B','C','D']
            k=keys[j//2]; re=(j%2==0)
            vals[k]=1.0 if re else 1j
            vals['Abar']=np.conj(vals['A']); vals['Bbar']=np.conj(vals['B'])
            vals['Cbar']=np.conj(vals['C']); vals['Dbar']=np.conj(vals['D'])
            F=sum(cF[n]*vals[n] for n in names)
            T=sum(cT[n]*vals[n] for n in names)
            syy=np.real(F); sxy=-np.imag(F); sxx=np.real(T)-syy
            M[:,j]=[sxx,syy,sxy]
        rows.append(M); rhs.append(target)
    return np.vstack(rows), np.concatenate(rhs)

# sample points in upper half plane, away from origin/crack
pts=[]
for x in [-2,-1.2,-0.5,0.3,0.8,1.5,2.5]:
    for y in [0.15,0.3,0.55,0.9]:
        pts.append((x,y))
Amat,bvec=build_system(c1,pts)
sol,res,rank,sv=np.linalg.lstsq(Amat,bvec,rcond=None)
print("lstsq residual",np.linalg.norm(Amat@sol-bvec)/np.linalg.norm(bvec),"rank",rank,"sv",sv)
A=sol[0]+1j*sol[1]; B=sol[2]+1j*sol[3]; C=sol[4]+1j*sol[5]; D=sol[6]+1j*sol[7]
print("A,B,C,D =",A,B,C,D)
# compare to expected from a+\bar b etc
a,b,c,d=c1
print("a+conj(b)=",a+np.conj(b)," c+conj(d)=",c+np.conj(d))

# evaluate stresses and displacement from fitted potentials
def eval_pot(A,B,C,D,x,y):
    z=x+1j*y; zc=np.conj(z)
    l=lam; lb=np.conj(l)
    phi=A*z**l+B*z**lb
    Phi=l*A*z**(l-1)+lb*B*z**(lb-1)
    Phip=l*(l-1)*A*z**(l-2)+lb*(lb-1)*B*z**(lb-2)
    Phic=np.conj(Phi); Phipc=np.conj(Phip)
    psi=(l+1)*C*z**l+(lb+1)*D*z**lb
    Psi=(l+1)*l*C*z**(l-1)+(lb+1)*lb*D*z**(lb-1)
    F=Phi+Phic+z*Phipc+np.conj(Psi)
    T=2*(Phi+Phic)
    syy=np.real(F); sxy=-np.imag(F); sxx=np.real(T)-syy
    E,nup,mu=Estar1,nuprime,mu1
    Dp=(kappa*phi-z*np.conj(Phi)-np.conj(psi))/(2*mu)
    return sxx,syy,sxy,np.real(Dp),np.imag(Dp)

print("stress check at sample points:")
for (x,y) in [(0.5,0.3),(1.0,0.2),(2.0,1.0),(-0.7,0.4)]:
    r=np.hypot(x,y); th=np.arctan2(y,x)
    f=c1@basis(th,s); fp=c1@dbasis(th,s); fpp=c1@ddbasis(th,s)
    Prr=(s+2)*f+fpp; Ptt=(s+2)*(s+1)*f; Prt=-(s+1)*fp; rs=r**s
    Sxx=np.real(rs*(Prr*np.cos(th)**2+Ptt*np.sin(th)**2-2*Prt*np.sin(th)*np.cos(th)))
    Syy=np.real(rs*(Prr*np.sin(th)**2+Ptt*np.cos(th)**2+2*Prt*np.sin(th)*np.cos(th)))
    Sxy=np.real(rs*((Prr-Ptt)*np.sin(th)*np.cos(th)+Prt*(np.cos(th)**2-np.sin(th)**2)))
    pp=eval_pot(A,B,C,D,x,y)
    print(" ",x,y,"airy",np.round([Sxx,Syy,Sxy],6),"pot",np.round(pp[:3],6),"diff",np.round(np.array(pp[:3])-np.array([Sxx,Syy,Sxy]),8))
