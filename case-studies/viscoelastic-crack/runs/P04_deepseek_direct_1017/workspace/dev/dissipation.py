import numpy as np

def fields_T(r, V, lam=0.1):
    """Compute E1,E2,T1,T2 on the upper-half interior grid (j=1..Ny-1, i=1..Nx-2)
    from solver output r.  U has shape (Ny-1, Nx-2)."""
    X = r['X']; U = r['U']; dY = r['dY']; Ny = U.shape[0]+1; Nx = len(X)
    # X derivative of U at interior i=1..Nx-2 (nonuniform central, one-sided at edges)
    E1 = np.zeros_like(U); E2 = np.zeros_like(U)
    for k in range(U.shape[1]):
        i = k+1
        if i-1 >= 1 and i+1 <= Nx-2:
            hm = X[i]-X[i-1]; hp = X[i+1]-X[i]
            E1[:,k] = np.pi*(-hp/(hm*(hm+hp))*U[:,k-1] + (hp-hm)/(hm*hp)*U[:,k] + hm/(hp*(hm+hp))*U[:,k+1])
        elif i-1 < 1:
            E1[:,k] = np.pi*(U[:,k+1]-U[:,k])/(X[i+1]-X[i])
        else:
            E1[:,k] = np.pi*(U[:,k]-U[:,k-1])/(X[i]-X[i-1])
    # Y derivative: j=1..Ny-1; top boundary U=a at Y=pi (need a from far field T2? pass separately)
    # We'll handle top with the value from r['far_a'] if present, else extrapolate.
    for j in range(U.shape[0]):
        if j+1 <= U.shape[0]-1:
            E2[j,:] = np.pi*(U[j+1,:]-U[j,:])/dY
        else:
            # top: use boundary value a from r
            a = r.get('a', None)
            if a is None:
                # extrapolate from last two
                E2[j,:] = np.pi*(U[j,:]-U[j-1,:])/dY
            else:
                E2[j,:] = np.pi*(a-U[j,:])/dY
    # T from ODE backward in X at each Y
    T1 = np.zeros_like(U); T2 = np.zeros_like(U)
    # far field at X=L (i=Nx-1): E1=0,E2=a,T=E
    a = r.get('a', 0.0)
    for j in range(U.shape[0]):
        for comp,E,T in ((0,E1,T1),(1,E2,T2)):
            e = E[j,:]
            t = np.zeros_like(e)
            # boundary value at i=Nx-1: far field
            # For comp 0 E1=0; comp 1 E2=a
            e_right = 0.0 if comp==0 else a
            t_right = e_right
            t_prev = t_right; e_prev = e_right
            for k in range(len(e)-1, -1, -1):
                i = k+1
                h = X[i+1]-X[i]
                R = e[k] - V*(e_prev-e[k])/h
                t[k] = (R + lam*V*t_prev/h)/(1.0+lam*V/h)
                t_prev = t[k]; e_prev = e[k]
            T[j,:] = t
    return E1,E2,T1,T2

def gamma_D_ratio(r, V, lam=0.1, delta_bar=0.2):
    if V == 0: return 0.0
    E1,E2,T1,T2 = fields_T(r,V,lam)
    X = r['X']; dY = r['dY']
    Xi = X[1:-1]
    integ = (T1-E1)**2 + (T2-E2)**2
    # integrate over X nonuniform (trapezoid) and Y uniform, upper half only
    Ix = np.zeros(integ.shape[0])
    for j in range(integ.shape[0]):
        Ix[j] = np.trapezoid(integ[j,:], Xi)
    I = np.trapezoid(Ix, np.arange(1,integ.shape[0]+1)*dY)
    pref = 2*lam/(np.pi*V*(1-lam)*delta_bar)
    return pref*I
