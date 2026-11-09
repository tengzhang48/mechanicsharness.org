"""
Reconsider the eigenvalue counting.

J-integral scaling: J = int (W cos - t . du/dx) r dtheta.
With sigma ~ r^{-p}, eps ~ sigma^n ~ r^{-pn}, W ~ sigma*eps ~ r^{-p(1+n)}.
u ~ eps*r ~ r^{-pn-p+1} = r^{1-p(n+1)}.
du/dx ~ u/r ~ r^{-p(n+1)}.
int(...) r dtheta: W cos * r ~ r^{-p(n+1)};  t.du/dx * r ~ r^{-p}*r^{-p(n+1)}*r = r^{-p(n+1)}.
So J ~ r^{-p(n+1)} * const(theta).  For J to be finite & r-independent: p(n+1)=1 -> p=1/(n+1).
"""
print("Dimensional argument: J = const(theta) * r^{-p(n+1)}, finite => p = 1/(n+1).")
print("So p is NOT an eigenvalue; p = 1/(n+1) exactly.")
print()
print("REVISED scalings (n>1 total-deformation) [p=1/(n+1)]:")
print(" sigma_e = sigma0 K r^{-p} sebar ;  eps_e = alpha eps0 K^n r^{-pn} sebar^n")
print(" pn = n/(n+1) = 1-p  =>  eps ~ r^{-(1-p)}")
print(" displacement u ~ eps*r ~ r^{1-(1-p)} = r^{p}")
print(" W = n/(n+1) sigma:eps ~ r^{-p}*r^{-(1-p)} = r^{-1}")
print(" => U(R) = int_0^R int W r dth dr ~ int_0^R (r^{-1}*r) dr = R")
print(" => U(R) is PROPORTIONAL TO R.  C = U(R)/(J R) is R-INDEPENDENT. Good.")
print()
print("Check J: J = int(W cos - t_i du_i/dx) r dth ; t~r^{-p}, du/dx~u/r~r^{p-1}")
print(" => t.du/dx ~ r^{-1}; W~r^{-1}; times r => r^0. J is R-independent. Good.")
print("So with CORRECT scalings everything is consistent: J is the constant,")
print("U(R) ~ R, C = U/(JR) is a pure number.")