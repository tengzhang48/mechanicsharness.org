import sympy as sp, pickle
k, s, G, c = sp.symbols('k s G c', positive=True)
nu, nuu, al = sp.symbols('nu nu_u alpha', positive=True)
m = sp.symbols('m', positive=True)
d = pickle.load(open('<opus-workspace>/code/symbols.pkl', 'rb'))
rep = {c: 1}
A11 = d['A11'].subs(rep)
A11m = sp.simplify(A11.subs(sp.sqrt(k**2 + s), m).subs(s, m**2 - k**2))
print('A11(m) =', sp.factor(A11m))
# guess: A11 = -G k/(1-nu_u) * [1 + X]
X = sp.simplify(sp.factor(A11m/(-G*k/(1-nuu)) - 1))
print('X =', sp.factor(X))
