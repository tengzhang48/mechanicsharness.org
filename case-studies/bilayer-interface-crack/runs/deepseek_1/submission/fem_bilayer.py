"""
Plane-strain high-order FEM for a bilayer strip with a semi-infinite interface crack.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ---------- 1D Lagrange basis on [-1,1], equispaced nodes ----------
def lagrange_1d(p, xi):
    nodes = np.linspace(-1.0, 1.0, p+1)
    xi = np.asarray(xi, dtype=float)
    n = xi.size
    L = np.zeros((p+1, n))
    dL = np.zeros((p+1, n))
    for i in range(p+1):
        num = np.ones(n)
        den = 1.0
        for j in range(p+1):
            if j == i:
                continue
            num *= (xi - nodes[j])
            den *= (nodes[i] - nodes[j])
        L[i] = num/den
        dsum = np.zeros(n)
        for k in range(p+1):
            if k == i:
                continue
            prod = np.ones(n)
            for j in range(p+1):
                if j == i or j == k:
                    continue
                prod *= (xi - nodes[j])
            dsum += prod
        dL[i] = dsum/den
    return L, dL


def graded_nodes(a, b, n, beta, cluster_at='a'):
    """n macro elements from a to b; geometric-like clustering at 'a' or 'b'."""
    t = np.linspace(0.0, 1.0, n+1)
    if cluster_at == 'a':
        s = t**beta
        return a + (b-a)*s
    elif cluster_at == 'b':
        s = 1.0 - (1.0-t)**beta
        return a + (b-a)*s
    else:
        raise ValueError


class Mesh:
    def __init__(self, h=1.0, H=2.0, L=20.0, p=3,
                 nxL=60, nxR=60, beta_x=2.5,
                 ny1=18, ny2=22, gamma_y=2.0,
                 crack=True):
        self.h, self.H, self.L, self.p = h, H, L, p
        self.crack = crack
        # macro x nodes
        xl = -L*((nxL-np.arange(nxL+1))/nxL)**beta_x  # -L ... 0
        xr = L*(np.arange(nxR+1)/nxR)**beta_x          # 0 ... L
        Xm = np.concatenate([xl[:-1], xr])             # nxL+nxR+1 macro nodes
        self.Xm = Xm
        self.nx = nxL + nxR
        self.i0 = nxL  # index of x=0 macro node
        # macro y nodes
        Y1m = h*(np.arange(ny1+1)/ny1)**gamma_y
        Y2m = -H*((ny2-np.arange(ny2+1))/ny2)**gamma_y
        self.Y1m, self.Y2m = Y1m, Y2m
        self.ny1, self.ny2 = ny1, ny2
        # fine 1D nodes
        self.x = self._fine(Xm, p)
        self.y1 = self._fine(Y1m, p)
        self.y2 = self._fine(Y2m, p)
        self.Nx = self.x.size - 1
        self.Ny1 = self.y1.size - 1
        self.Ny2 = self.y2.size - 1
        assert abs(self.x[self.i0*p]) < 1e-14
        # global nodes
        coords = []
        node1 = -np.ones((self.Nx+1, self.Ny1+1), dtype=int)
        node2 = -np.ones((self.Nx+1, self.Ny2+1), dtype=int)
        for i in range(self.Nx+1):
            for j in range(self.Ny1+1):
                node1[i,j] = len(coords)
                coords.append((self.x[i], self.y1[j]))
        for i in range(self.Nx+1):
            for j in range(self.Ny2+1):
                if j == self.Ny2 and (self.x[i] > -1e-12 or not crack):
                    # interface: merge with layer-1 y=0 node (bonded, or fully bonded)
                    node2[i,j] = node1[i,0]
                else:
                    node2[i,j] = len(coords)
                    coords.append((self.x[i], self.y2[j]))
        # when not cracked, x=0 interface also merged (condition above includes x>=0? x>1e-14 excludes exactly 0)
        self.node1, self.node2 = node1, node2
        self.coords = np.array(coords, dtype=float)
        self.nn = self.coords.shape[0]
        self.ndof = 2*self.nn
        # element connectivity lists
        self.elems = []  # (layer, e, f, global_ids[p+1,p+1])
        for e in range(self.nx):
            for f in range(self.ny1):
                ids = node1[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
                self.elems.append((1,e,f,ids))
            for f in range(self.ny2):
                ids = node2[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
                self.elems.append((2,e,f,ids))
        self.nelem = len(self.elems)

    @staticmethod
    def _fine(macro, p):
        out = [macro[0]]
        for k in range(len(macro)-1):
            a, b = macro[k], macro[k+1]
            out.extend(list(a + (b-a)*np.arange(1,p+1)/p))
        return np.array(out)


def D_plane_strain(E, nu):
    mu = E/(2*(1+nu))
    lam = E*nu/((1+nu)*(1-2*nu))
    return np.array([[lam+2*mu, lam, 0.0],
                     [lam, lam+2*mu, 0.0],
                     [0.0, 0.0, mu]])


class Solver:
    def __init__(self, mesh, E1=9.0, E2=1.0, nu1=1/3, nu2=1/3):
        self.m = mesh
        self.D1 = D_plane_strain(E1, nu1)
        self.D2 = D_plane_strain(E2, nu2)
        self.p = mesh.p
        # Gauss rule
        self.ng = self.p + 1
        self.gx, self.gw = np.polynomial.legendre.leggauss(self.ng)
        # shape functions at Gauss points (2D)
        Lg, dLg = lagrange_1d(self.p, self.gx)
        self.Lg, self.dLg = Lg, dLg

    def assemble(self):
        m = self.m
        p, ng = self.p, self.ng
        nloc = (p+1)**2
        ndof_e = 2*nloc
        rows = []
        cols = []
        vals = []
        # precompute 2D shape function values/derivatives at Gauss points
        Lg, dLg = self.Lg, self.dLg
        for (layer, e, f, ids) in m.elems:
            xa, xb = m.x[e*p], m.x[(e+1)*p]
            if layer == 1:
                ya, yb = m.y1[f*p], m.y1[(f+1)*p]
                D = self.D1
            else:
                ya, yb = m.y2[f*p], m.y2[(f+1)*p]
                D = self.D2
            dx, dy = xb-xa, yb-ya
            # local dof map
            ldof = np.empty(ndof_e, dtype=int)
            k = 0
            for b in range(p+1):
                for a in range(p+1):  # order: node (a,b)
                    gid = ids[a,b]
                    ldof[k] = 2*gid
                    ldof[k+1] = 2*gid+1
                    k += 2
            Ke = np.zeros((ndof_e, ndof_e))
            # note: ids[a,b] indexing: ids has shape (p+1,p+1) with first index x, second y
            # but in the loop above ids[a,b] uses a=x index, b=y index. local dof k order:
            # b outer, a inner => matches B construction below
            for gi in range(ng):
                xi = self.gx[gi]
                wxi = self.gw[gi]
                for gj in range(ng):
                    eta = self.gx[gj]
                    weta = self.gw[gj]
                    N = np.outer(Lg[:,gi], Lg[:,gj])  # N[a,b]? outer first index a, second b
                    dNdxi = np.outer(dLg[:,gi], Lg[:,gj])
                    dNdeta = np.outer(Lg[:,gi], dLg[:,gj])
                    dNdx = dNdxi * (2.0/dx)
                    dNdy = dNdeta * (2.0/dy)
                    B = np.zeros((3, ndof_e))
                    k = 0
                    for b in range(p+1):
                        for a in range(p+1):
                            B[0,k] = dNdx[a,b]
                            B[1,k+1] = dNdy[a,b]
                            B[2,k] = dNdy[a,b]
                            B[2,k+1] = dNdx[a,b]
                            k += 2
                    Ke += (wxi*weta)*(dx/2.0)*(dy/2.0)*(B.T @ D @ B)
            rr, cc = np.meshgrid(ldof, ldof, indexing='ij')
            rows.append(rr.ravel())
            cols.append(cc.ravel())
            vals.append(Ke.ravel())
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        vals = np.concatenate(vals)
        self.K = sp.coo_matrix((vals, (rows, cols)), shape=(m.ndof, m.ndof)).tocsc()
        return self.K

    def load_left(self, P, M):
        """Consistent nodal forces on x=-L from the prescribed beam tractions."""
        m = self.m
        p, ng = self.p, self.ng
        f = np.zeros(m.ndof)
        # Gauss rule for boundary (p+1 points)
        gx, gw = np.polynomial.legendre.leggauss(self.ng)
        Lg, dLg = lagrange_1d(p, gx)
        # layer 1 and layer 2 boundary element edges
        for layer, ny, yarr in [(1, m.ny1, m.y1), (2, m.ny2, m.y2)]:
            ids = (m.node1 if layer==1 else m.node2)[0, :]  # x index 0, all y
            for f_e in range(ny):
                ya, yb = yarr[f_e*p], yarr[(f_e+1)*p]
                dy = yb-ya
                yc = 0.5*(ya+yb)
                # nodal ids for this edge, b index
                edge_ids = ids[f_e*p:(f_e+1)*p+1]
                for gi in range(self.ng):
                    eta = gx[gi]; w = gw[gi]
                    y = yc + 0.5*dy*eta
                    if layer == 1:
                        sig11 = -P/m.h - 12.0*M*(y - m.h/2.0)/m.h**3
                    else:
                        M2 = M + P*(m.h+m.H)/2.0
                        sig11 = P/m.H + 12.0*M2*(y + m.H/2.0)/m.H**3
                    tx = -sig11  # external traction on outward normal (-x)
                    for b in range(p+1):
                        f[2*edge_ids[b]] += w * (dy/2.0) * Lg[b,gi] * tx
        return f

    def solve(self, P=0.0, M=0.0, fix=None):
        m = self.m
        if not hasattr(self, 'K'):
            self.assemble()
        f = self.load_left(P, M)
        # constraints: fix rigid body modes at right boundary
        # node at (L,-H) = node2[Nx,0]; node at (L,h)=node1[Nx,Ny1]
        nbr = m.node2[m.Nx, 0]
        ntr = m.node1[m.Nx, m.Ny1]
        fixed = [2*nbr, 2*nbr+1, 2*ntr]
        free = np.ones(m.ndof, dtype=bool)
        free[fixed] = False
        u = np.zeros(m.ndof)
        u[free] = spla.spsolve(self.K[free][:,free], f[free])
        return u

    # ---------- postprocessing ----------
    def disp(self, u, gid):
        return np.array([u[2*gid], u[2*gid+1]])

    def cod(self, u):
        """Return r>0, delta2+i delta1 on the crack faces (upper - lower)."""
        m = self.m
        i0 = m.i0*m.p
        rs = []
        ds = []
        for i in range(i0):
            r = -m.x[i]
            if r <= 0:
                continue
            gu = m.node1[i,0]
            gl = m.node2[i,m.Ny2]
            if gu == gl:
                continue
            du = self.disp(u,gu)
            dl = self.disp(u,gl)
            d = du-dl
            rs.append(r)
            ds.append(d[1] + 1j*d[0])  # delta2 + i delta1
        return np.array(rs), np.array(ds)

    def stress_interface(self, u, x, side='avg'):
        """Stress at interface y=0, x>0, from adjacent elements above/below."""
        m = self.m
        p = self.p
        # find macro/fine element index
        i = np.searchsorted(m.x, x, side='right') - 1
        if i < 0 or i >= m.Nx:
            raise ValueError('x outside mesh')
        e = i//p
        # local xi
        xa, xb = m.x[e*p], m.x[(e+1)*p]
        xi = 2.0*(x-xa)/(xb-xa) - 1.0
        out = []
        # above: layer1 f=0, eta=-1
        if side in ('avg','above'):
            sig = self._stress_elem(u, 1, e, 0, xi, -1.0)
            out.append(sig)
        if side in ('avg','below'):
            sig = self._stress_elem(u, 2, e, m.ny2-1, xi, 1.0)
            out.append(sig)
        return np.mean(out, axis=0)

    def _stress_elem(self, u, layer, e, f, xi, eta):
        m = self.m
        p = self.p
        if layer == 1:
            D = self.D1
            ids = m.node1[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
            ya, yb = m.y1[f*p], m.y1[(f+1)*p]
        else:
            D = self.D2
            ids = m.node2[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
            ya, yb = m.y2[f*p], m.y2[(f+1)*p]
        xa, xb = m.x[e*p], m.x[(e+1)*p]
        dx, dy = xb-xa, yb-ya
        L, dL = lagrange_1d(p, np.array([xi]))
        M, dM = lagrange_1d(p, np.array([eta]))
        L, dL, M, dM = L[:,0], dL[:,0], M[:,0], dM[:,0]
        dNdx = np.outer(dL, M)*(2.0/dx)
        dNdy = np.outer(L, dM)*(2.0/dy)
        ue = np.zeros(2*(p+1)**2)
        k=0
        for b in range(p+1):
            for a in range(p+1):
                gid=ids[a,b]
                ue[k]=u[2*gid]; ue[k+1]=u[2*gid+1]; k+=2
        eps = np.zeros(3)
        k=0
        for b in range(p+1):
            for a in range(p+1):
                eps[0] += dNdx[a,b]*ue[k]
                eps[1] += dNdy[a,b]*ue[k+1]
                eps[2] += dNdy[a,b]*ue[k] + dNdx[a,b]*ue[k+1]
                k+=2
        sig = D @ eps
        return sig  # [s11,s22,s12]

    def _grad_elem(self, u, layer, e, f, xi, eta):
        """Return (eps(3), grad_u (2x2)) at a local point."""
        m = self.m
        p = self.p
        if layer == 1:
            D = self.D1
            ids = m.node1[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
            ya, yb = m.y1[f*p], m.y1[(f+1)*p]
        else:
            D = self.D2
            ids = m.node2[e*p:(e+1)*p+1, f*p:(f+1)*p+1]
            ya, yb = m.y2[f*p], m.y2[(f+1)*p]
        xa, xb = m.x[e*p], m.x[(e+1)*p]
        dx, dy = xb-xa, yb-ya
        L, dL = lagrange_1d(p, np.array([xi]))
        M, dM = lagrange_1d(p, np.array([eta]))
        L, dL, M, dM = L[:,0], dL[:,0], M[:,0], dM[:,0]
        dNdx = np.outer(dL, M)*(2.0/dx)
        dNdy = np.outer(L, dM)*(2.0/dy)
        ue = np.zeros(2*(p+1)**2)
        k=0
        for b in range(p+1):
            for a in range(p+1):
                gid=ids[a,b]
                ue[k]=u[2*gid]; ue[k+1]=u[2*gid+1]; k+=2
        eps = np.zeros(3); g = np.zeros((2,2))
        k=0
        for b in range(p+1):
            for a in range(p+1):
                eps[0] += dNdx[a,b]*ue[k]
                eps[1] += dNdy[a,b]*ue[k+1]
                eps[2] += dNdy[a,b]*ue[k] + dNdx[a,b]*ue[k+1]
                g[0,0] += dNdx[a,b]*ue[k]
                g[0,1] += dNdy[a,b]*ue[k]
                g[1,0] += dNdx[a,b]*ue[k+1]
                g[1,1] += dNdy[a,b]*ue[k+1]
                k+=2
        sig = D @ eps
        return sig, eps, g

    def J_line(self, u, x0):
        """J-integral along vertical line x=x0<0 from y=-H to h (both layers)."""
        m = self.m
        # find element column containing x0
        i = np.searchsorted(m.x, x0, side='right') - 1
        if i < 0 or i >= m.Nx:
            raise ValueError('x0 outside mesh')
        e = i//self.p
        xa, xb = m.x[e*self.p], m.x[(e+1)*self.p]
        xi = 2.0*(x0-xa)/(xb-xa) - 1.0
        # Gauss in y
        gx, gw = np.polynomial.legendre.leggauss(self.ng)
        J = 0.0
        for layer, ny, yarr in [(2, m.ny2, m.y2), (1, m.ny1, m.y1)]:
            for f in range(ny):
                ya, yb = yarr[f*self.p], yarr[(f+1)*self.p]
                dy = yb-ya; yc=0.5*(ya+yb)
                for gi in range(self.ng):
                    eta = gx[gi]; w=gw[gi]
                    sig, eps, g = self._grad_elem(u, layer, e, f, xi, eta)
                    W = 0.5*(sig[0]*eps[0]+sig[1]*eps[1]+2*sig[2]*eps[2])
                    # n=(-1,0): integrand = -W + sig11*u1x + sig12*u2x
                    integ = -W + sig[0]*g[0,0] + sig[2]*g[1,0]
                    J += w*(dy/2.0)*integ
        return J

    def eval_point(self, u, x, y):
        """Return (sigma[3], eps[3], grad[2,2]) at an arbitrary point."""
        m = self.m
        p = self.p
        if y >= 0.0:
            layer=1; yarr=m.y1; ny=m.ny1; ids=m.node1
        else:
            layer=2; yarr=m.y2; ny=m.ny2; ids=m.node2
        i = np.searchsorted(m.x, x, side='right') - 1
        i = max(0, min(i, m.Nx-1))
        e = i//p
        xa, xb = m.x[e*p], m.x[(e+1)*p]
        xi = 2.0*(x-xa)/(xb-xa) - 1.0
        j = np.searchsorted(yarr, y, side='right') - 1
        j = max(0, min(j, ny-1))
        f = j//p
        ya, yb = yarr[f*p], yarr[(f+1)*p]
        eta = 2.0*(y-ya)/(yb-ya) - 1.0
        sig, eps, g = self._grad_elem(u, layer, e, f, xi, eta)
        return sig, eps, g

    def interaction_integral(self, u, aux, K1a=1.0, K2a=0.0, R=0.5, ntheta=120):
        """Contour interaction integral M for aux field with K=(K1a,K2a)."""
        gx, gw = np.polynomial.legendre.leggauss(ntheta)
        M = 0.0
        for lo,hi in [(-np.pi,0.0),(0.0,np.pi)]:
            th = 0.5*(hi-lo)*gx + 0.5*(hi+lo)
            w = 0.5*(hi-lo)*gw
            for t,wt in zip(th,w):
                x=R*np.cos(t); y=R*np.sin(t)
                sig,eps,g = self.eval_point(u,x,y)
                siga,ua,ga = aux.fields(x,y,K1a,K2a)
                n=np.array([np.cos(t),np.sin(t)])
                tvec=np.array([sig[0]*n[0]+sig[2]*n[1], sig[2]*n[0]+sig[1]*n[1]])
                tveca=np.array([siga[0]*n[0]+siga[2]*n[1], siga[2]*n[0]+siga[1]*n[1]])
                # W_int = sigma_ij eps^aux_ij
                epsa = np.array([ga[0,0], ga[1,1], 0.5*(ga[0,1]+ga[1,0])])
                Wint = sig[0]*epsa[0] + sig[1]*epsa[1] + 2.0*sig[2]*epsa[2]
                integ = Wint*n[0] - (tvec[0]*ga[0,0]+tvec[1]*ga[1,0]) \
                                 - (tveca[0]*g[0,0]+tveca[1]*g[1,0])
                M += wt*R*integ
        return M

    def interaction_line(self, u, aux, K1a=1.0, K2a=0.0, x0=-5.0):
        """Interaction integral on a vertical line x=x0<0, from y=-H to h."""
        m = self.m
        i = np.searchsorted(m.x, x0, side='right') - 1
        if i < 0 or i >= m.Nx:
            raise ValueError('x0 outside mesh')
        e = i//self.p
        xa, xb = m.x[e*self.p], m.x[(e+1)*self.p]
        xi = 2.0*(x0-xa)/(xb-xa) - 1.0
        gx, gw = np.polynomial.legendre.leggauss(self.ng)
        M = 0.0
        for layer, ny, yarr in [(2, m.ny2, m.y2), (1, m.ny1, m.y1)]:
            for f in range(ny):
                ya, yb = yarr[f*self.p], yarr[(f+1)*self.p]
                dy = yb-ya; yc=0.5*(ya+yb)
                for gi in range(self.ng):
                    eta = gx[gi]; w=gw[gi]; y=yc+0.5*dy*eta
                    sig, eps, g = self._grad_elem(u, layer, e, f, xi, eta)
                    siga, ua, ga = aux.fields(x0, y, K1a, K2a)
                    n=np.array([-1.0,0.0])
                    tvec=np.array([sig[0]*n[0]+sig[2]*n[1], sig[2]*n[0]+sig[1]*n[1]])
                    tveca=np.array([siga[0]*n[0]+siga[2]*n[1], siga[2]*n[0]+siga[1]*n[1]])
                    epsa=np.array([ga[0,0], ga[1,1], 0.5*(ga[0,1]+ga[1,0])])
                    Wint=sig[0]*epsa[0]+sig[1]*epsa[1]+2.0*sig[2]*epsa[2]
                    integ=Wint*n[0]-(tvec[0]*ga[0,0]+tvec[1]*ga[1,0])-(tveca[0]*g[0,0]+tveca[1]*g[1,0])
                    M += w*(dy/2.0)*integ
        return M

    def interaction_domain(self, u, aux, K1a=1.0, K2a=0.0, R1=0.1, R2=0.5):
        """Domain form of the interaction integral with a radial q-function."""
        m = self.m
        p, ng = self.p, self.ng
        gx, gw = np.polynomial.legendre.leggauss(ng)
        Lg, dLg = lagrange_1d(p, gx)
        M = 0.0
        for layer, e, f, ids in m.elems:
            if layer==1:
                xa,xb = m.x[e*p], m.x[(e+1)*p]
                ya,yb = m.y1[f*p], m.y1[(f+1)*p]
            else:
                xa,xb = m.x[e*p], m.x[(e+1)*p]
                ya,yb = m.y2[f*p], m.y2[(f+1)*p]
            # quick distance check
            ptsx=[xa,xb]; ptsy=[ya,yb]
            dmin=1e9; dmax=0
            for xx in ptsx:
                for yy in ptsy:
                    rr=np.hypot(xx,yy); dmin=min(dmin,rr); dmax=max(dmax,rr)
            if dmax < R1 or dmin > R2:
                continue
            dx,dy=xb-xa,yb-ya
            xc,yc=0.5*(xa+xb),0.5*(ya+yb)
            for gi in range(ng):
                for gj in range(ng):
                    xi=gx[gi]; eta=gx[gj]; w=gw[gi]*gw[gj]
                    x=xc+0.5*dx*xi; y=yc+0.5*dy*eta
                    r=np.hypot(x,y)
                    if r < R1 or r > R2:
                        continue
                    q=(R2-r)/(R2-R1); dqdr=-1.0/(R2-R1)
                    qx=dqdr*x/r; qy=dqdr*y/r
                    sig,eps,g = self._grad_elem(u, layer, e, f, xi, eta)
                    siga,ua,ga = aux.fields(x,y,K1a,K2a)
                    epsa=np.array([ga[0,0],ga[1,1],0.5*(ga[0,1]+ga[1,0])])
                    Wint=sig[0]*epsa[0]+sig[1]*epsa[1]+2.0*sig[2]*epsa[2]
                    T1=(sig[0]*ga[0,0]+sig[2]*ga[1,0])*qx + (sig[2]*ga[0,0]+sig[1]*ga[1,0])*qy
                    T2=(siga[0]*g[0,0]+siga[2]*g[1,0])*qx + (siga[2]*g[0,0]+siga[1]*g[1,0])*qy
                    integ=T1+T2-Wint*qx
                    M += w*(dx/2.0)*(dy/2.0)*integ
        return M
