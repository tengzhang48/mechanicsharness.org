import numpy as np


def graded_nodes(N, p=3.0, xc=0.0):
    """Graded nodes on [0,1] clustering near x=0 using the mapping
       x(eta) = 1 - (1-eta)^?  ... we use x = eta^p with eta uniform in [0,1].
       Returns the FULL node array of length N+2 (endpoints 0 and 1)."""
    eta = np.linspace(0.0, 1.0, N + 2)
    x = eta ** p
    x[0] = 0.0
    x[-1] = 1.0
    return x


def graded_nodes_geo(N, r=1.0, Lcluster=0.5):
    """Geometric grading near x=0 then uniform.  First n0 points geometric
    with ratio r (spacing), covering [0, Lcluster], rest uniform to 1."""
    # split N+2 nodes: n0 in cluster, rest uniform
    n0 = (N + 2) // 2
    rest = (N + 2) - n0
    # geometric spacing d_k = d0 * r^k ; choose d0 so sum = Lcluster
    if abs(r - 1.0) < 1e-12:
        d = np.full(n0, Lcluster / n0)
    else:
        k = np.arange(n0)
        d0 = Lcluster * (r - 1) / (r ** n0 - 1)
        d = d0 * r ** k
    xc = np.concatenate(([0.0], np.cumsum(d)))
    xc = xc / xc[-1] * Lcluster
    xr = np.linspace(Lcluster, 1.0, rest)
    x = np.concatenate((xc, xr[1:]))
    x[0] = 0.0
    x[-1] = 1.0
    return x


if __name__ == "__main__":
    for p in [1, 2, 3, 4]:
        x = graded_nodes(20, p)
        print("p=%.0f min spacing %.4e  first nodes %s" % (p, np.diff(x).min(), np.round(x[:6], 4)))