"""Independent interface-traction extraction from a solved FE model."""
import numpy as np
import bilayer_fe as b


def interface_stress_K(res):
    x, nx, nb = res['x'], res['nx'], res['nb']
    tip = nx // 2
    cells = np.arange(tip, nx - 1)
    xc = 0.5 * (x[cells] + x[cells + 1])
    Ks = np.zeros((len(cells), 2), complex)
    U = res['U']
    for material in ('bottom', 'top'):
        E = b.E2 if material == 'bottom' else b.E1
        D = E / ((1 + b.NU) * (1 - 2 * b.NU)) * np.array(
            [[1 - b.NU, b.NU, 0], [b.NU, 1 - b.NU, 0], [0, 0, (1 - 2*b.NU)/2]])
        eta0 = 1.0 if material == 'bottom' else -1.0
        vals = np.zeros((len(cells), 2), complex)
        for jj, ix in enumerate(cells):
            if material == 'bottom':
                n0 = (nb - 2) * nx + ix
                conn = np.array([n0, n0+1, n0+1+nx, n0+nx])
                hy = res['yb'][-1] - res['yb'][-2]
            else:
                n0 = ix
                conn = res['top_map'][np.array([n0,n0+1,n0+1+nx,n0+nx])]
                hy = res['yt'][1] - res['yt'][0]
            hx = x[ix+1]-x[ix]
            dofs = np.empty(8,dtype=int);dofs[0::2]=2*conn;dofs[1::2]=2*conn+1
            # sigma_yy and sigma_xy vary linearly in xi for this Q4 row at
            # the interface midpoint, so use the center value xi=0.
            xi,eta=0.0,eta0
            dN_dxi=.25*np.array([-(1-eta),(1-eta),(1+eta),-(1+eta)])
            dN_deta=.25*np.array([-(1-xi),-(1+xi),(1+xi),(1-xi)])
            dx=dN_dxi*2/hx;dy=dN_deta*2/hy
            B=np.zeros((3,8));B[0,0::2]=dx;B[1,1::2]=dy;B[2,0::2]=dy;B[2,1::2]=dx
            stress=D@B@U[dofs]
            vals[jj]=stress[1]+1j*stress[2]
        Ks += np.sqrt(2*np.pi*xc[:,None])*vals*np.exp(-1j*b.EPS*np.log(xc[:,None]))/2
    return xc,Ks


def fits(r,K):
    for rmax in [.02,.04,.08,.16,.32,.64]:
        mask=r<rmax
        z=np.sqrt(r[mask]/rmax)
        for deg in [1,2,3]:
            c=np.linalg.lstsq(np.vander(z,deg+1,increasing=True),K[mask],rcond=None)[0]
            print('traction',rmax,deg,[(abs(v),np.degrees(np.angle(v)),v.real,v.imag) for v in c[0]])


if __name__=='__main__':
    res=b.solve_model(first=3e-5,growth=1.10,largest_x=.12,largest_y=.10)
    r,K=interface_stress_K(res)
    fits(r,K)
    np.savez('traction_medium.npz',r=r,K=K)
