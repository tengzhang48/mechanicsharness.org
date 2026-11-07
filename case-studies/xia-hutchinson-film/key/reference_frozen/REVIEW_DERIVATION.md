# Independent derivation (checker)

Authoritative energy (PROBLEM.md, not a guessed reading of XH2000 Eq. 6):

    Delta_Pi = int [ h sigma0 eps_aa + (h/2) Delta_sigma:eps + (k/2) |u|^2 ] dA

Ebar = E/(1-nu^2). Plane-stress constitutive law as stated. Total traction sigma n = 0 on faces.

## Strip equilibrium

Far from tips, x1-independent. Mode-I: u1=0, u2=u(x2). Then eps22=u', Delta_sigma22=Ebar u',
Delta_sigma11=nu Ebar u', Delta_sigma12=0. Equilibrium h Ebar u'' = k u gives

    ell = sqrt(Ebar h / k),   u'' = u/ell^2.

Total-traction faces: sigma22=0 => u' = -sigma0/Ebar at every crack face (independent of normal sign).
Midline of a spacing-S cell: u2=0 by reflection, sigma12=0.

Dimensionless v = Ebar u / (sigma0 ell), eta=x2/ell, a=S/(2 ell):
v''=v on a full cell (0, 2a) with natural conditions v'(0)=v'(2a)=-1.
Exact: v(eta)=sinh(a-eta)/cosh(a) on (0,a) and odd about eta=a.

## Translating-tip energy (steady)

psi dx2 = ell (h sigma0^2 / Ebar) [ v' + (v'^2 + v^2)/2 ] d eta
I(S) = int_0^{S/ell} [...] d eta = -tanh(S/(2 ell))  (first integral (v')^2 - v^2 = sech^2(a))

Case A: period H, one moving tip, ahead u=0 (I=0), behind spacing-H wake.
G_A h = -I(H)*ell*(h sigma0^2/Ebar)  =>  gamma_A = G_A Ebar/(sigma0^2 ell) = tanh(H/(2 ell))

Case B: period 2H, one new moving tip; ahead S=2H, behind two S=H cells.
G_B h = [I(2H) wait: Pi_ahead - Pi_behind]
     = ell (h sigma0^2/Ebar) * ( -I_behind_over_2H + I_ahead_over_2H wait)
From Pi_bar(S)*width: G_B h = (Pi_bar(2H) - Pi_bar(H))*2H
gamma_B = 2 tanh(H/(2 ell)) - tanh(H/ell) = 2 G'(H) - G'(2H)

KII=0 by reflection. G=KI^2/E => kI = KI/(sigma0 sqrt(ell)) = sqrt((1-nu^2) gamma)

Limits at fixed h,E,nu,k (fixed ell):
H/ell->0: gamma_A ~ H/(2ell) (1 - (H/ell)^2/12), gamma_B ~ (H/ell)^3/4
H/ell->inf: gamma_A, gamma_B -> 1; sigma_c -> sqrt(Gc Ebar / ell)

Source: XH2000 Eqs. 11-12 calibrate k so isolated G matches Beuth (1); Eq. 27 is the isolated kI.
Eq. 6 typesetting is ambiguous; PROBLEM.md functional is used. Eqs. 28/30 OCR is garbled;
the paper's identity G=2G'(H)-G'(2H) plus Eq. 28 tanh form recover the target formulas.
This is a cosmetic/source-transcription issue, not a target-quantity error.
