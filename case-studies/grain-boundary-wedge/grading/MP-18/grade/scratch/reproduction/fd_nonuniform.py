import numpy as np


def build_D1D2(nodes):
    """Nonuniform D1, D2 on interior nodes nodes[1..M], full grid nodes[0]=0,
    nodes[M+1]=1.  BCs: s(0)=0 (left value 0), s'(1)=s''(1)=0 (right).
    Stencils mirror march.py's uniform D1/D2, generalised to variable spacing:

      interior node i (1<=i<=M-2):
        D1[i] = [-hp/(hm(hm+hp)), (hp-hm)/(hm hp), hm/(hp(hm+hp))]
        D2[i] = [ 2/(hm(hm+hp)),  -2/(hm hp),      2/(hp(hm+hp)) ]

      first interior node (i=0, left neighbour is x=0 with value 0):
        D1[0] = [ 1/hm ,  -hm/(hp(hm+hp)) ]   (2nd-order one-sided, s0=0)
        D2[0] = [ -2/(hm hp),  2/(hp(hm+hp)) ] (s0=0)

      last interior node (i=M-1, right neighbour is x=1):
        s'(1)=0 -> ghost s_{M+1} = s_{M-1} (reflection about x=1 for odd BC),
        s''(1)=0 -> 2nd-difference at x=1 is zero -> ghost for D2.
        We implement the two right-end BCs as the standard march.py pattern:
          D1[M-1] = -1/(2 hm), +1/(2 hm)   (central with ghost = left neighbour)
          D2[M-1] =  2/hm^2, -2/hm^2
        which is exactly march.py's D1(N-1),D2(N-1) generalised.
    """
    x = np.asarray(nodes, float)
    M = len(x) - 2
    h = np.diff(x)                      # h[k]=x[k+1]-x[k], k=0..M
    D1 = np.zeros((M, M))
    D2 = np.zeros((M, M))
    # interior rows i=1..M-2
    for i in range(1, M - 1):
        fi = i + 1
        hm = h[fi - 1]
        hp = h[fi]
        D1[i, i - 1] = -hp / (hm * (hm + hp))
        D1[i, i] = (hp - hm) / (hm * hp)
        D1[i, i + 1] = hm / (hp * (hm + hp))
        D2[i, i - 1] = 2.0 / (hm * (hm + hp))
        D2[i, i] = -2.0 / (hm * hp)
        D2[i, i + 1] = 2.0 / (hp * (hm + hp))
    # first row (i=0): mirror march.D1 (s_0=0 used as ghost, centred)
    hm = h[0]
    hp = h[1]
    D1[0, 0] = 0.0
    if M > 1:
        D1[0, 1] = 1.0 / (2 * hm)      # (s_2 - s_0)/(2 hm), s_0 = 0
    D2[0, 0] = -2.0 / (hm * hp)
    if M > 1:
        D2[0, 1] = 2.0 / (hp * (hm + hp))
    # last row (i=M-1): mirror march.D1/D2 (ghost beyond x=1 via BC)
    hm = h[M - 1]
    D1[M - 1, M - 2] = -1.0 / (2 * hm)
    D1[M - 1, M - 1] = 0.0
    D2[M - 1, M - 2] = 2.0 / hm ** 2
    D2[M - 1, M - 1] = -2.0 / hm ** 2
    return D1, D2


if __name__ == "__main__":
    N = 30
    xn = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], xn, [1.0]))
    D1, D2 = build_D1D2(nodes)
    # compare with march on uniform
    import march
    print("D1 match march:", np.allclose(D1, march.D1(N)))
    print("D2 match march:", np.allclose(D2, march.D2(N)))
    f = xn ** 3
    print("D1 f vs 3x^2 max err:", np.max(np.abs(D1 @ f - 3 * xn ** 2)))
    print("D2 f vs 6x  max err:", np.max(np.abs(D2 @ f - 6 * xn)))
    g = xn ** 2 - 2 * xn
    print("D1 g last (s'(1)=0):", (D1 @ g)[-1], " D2 g last (s''=2):", (D2 @ g)[-1])