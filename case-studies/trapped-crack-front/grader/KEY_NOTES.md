# Key notes for the grader (never given to a solver)

The reference values themselves are in `score.md` (column "Reference").

## 1. The right first-order theory

Source: H. Gao and J. R. Rice (1986), J. Appl. Mech. 53, 774-778, Eq. 15; mode 1 from J. R. Rice (1985). For a
half-plane crack whose front is x = a(z), with straight-front factors K1^0, K2^0, K3^0 that are uniform and do not
depend on the position of the front, and with the frame of the task (n normal to the front, t = n x e_y):

    K1(z) = K1^0 + p1 K1^0 l(z)
    K2(z) = K2^0 + p2 K2^0 l(z) - b2 K3^0 a'(z)
    K3(z) = K3^0 + p3 K3^0 l(z) + b3 K2^0 a'(z)

    l(z) = (1/2 pi) PV integral of (a(z') - a(z))/(z' - z)^2 dz'        (on cos(kz): l = -(|k|/2) cos(kz))

| Coefficient | Meaning | Right expression | Value at nu = 0.3 |
|---|---|---|---|
| p1 | nonlocal weight, mode 1 | 1 | 1 |
| p2 | nonlocal weight, mode 2 | (2 - 3 nu)/(2 - nu) | 11/17 = 0.647059 |
| p3 | nonlocal weight, mode 3 | (2 + nu)/(2 - nu) | 23/17 = 1.352941 |
| b2 | slope coupling into mode 2 (sign: minus K3^0 a') | 2/(2 - nu) | 20/17 = 1.176471 |
| b3 | slope coupling into mode 3 (sign: plus K2^0 a') | 2 (1 - nu)/(2 - nu) | 14/17 = 0.823529 |
| m3 | weight of K3^2 in G, relative to K1^2 (Irwin) | 1/(1 - nu) | 1.428571 |

Energy release rate: G = (1 - nu^2)/E (K1^2 + K2^2 + m3 K3^2). With the right coefficients the slope terms cancel in
G to first order (m3 b3 = b2), so G(z) = G^0 [1 + (2/rho) l(z)] with

    rho = (1 + r2^2 + m3 r3^2)/(p1 + p2 r2^2 + m3 p3 r3^2) = 0.920795      (r2 = 0.6, r3 = 0.8).

## 2. The right use of it

- Part A, a = abar + A cos(2 pi z/lambda): l = -pi (A/lambda) cos, a' = -2 pi (A/lambda) sin. So C1 = -pi p1,
  C2 = -pi p2 r2, S2 = +2 pi b2 r3, C3 = -pi p3 r3, S3 = -2 pi b3 r2, CG = -2 pi/rho, S1 = SG = 0.
- Part B: the mean of G = Gamma gives G^0 = Gamma_0, so K10_scaled = 1/sqrt(1 + r2^2 + m3 r3^2). With
  g(s) = 2 sum q^n cos(2 pi n s), the equilibrium equation is l = eps (rho/2) g, and
  f(s) = (rho/2 pi) ln(1 - 2 q cos 2 pi s + q^2). The front lags in the tough bands: f(0) < 0, f(1/2) > 0.
  f'(s) = rho S(s), S = 2 q sin(2 pi s)/(1 - 2 q cos 2 pi s + q^2);
  kappa1 = p1 (rho/2) g, kappa2 = p2 r2 (rho/2) g - b2 r3 rho S, kappa3 = p3 r3 (rho/2) g + b3 r2 rho S.

## 3. Usual mistakes

| Mistake | Kind | How it shows |
|---|---|---|
| Slope couplings 1 and 1 ("rotate the shear tractions into the tilted frame"), often with the mode 1 weight for all modes | theory | S2, S3 and the shear tables wrong; the slope terms no longer cancel in G, so the front gets an odd part |
| No nonlocal term at all (Meade and Keer 1984, quoted as Eq. 16 in Gao and Rice) | theory | C1, C2, C3 zero |
| No slope coupling | theory | S2 = S3 = 0; kappa2, kappa3 even in s |
| Mode 3 with weight 1 in G | theory | K10_scaled, CG and the front wrong |
| Nonlocal term doubled (1/pi for 1/(2 pi)) | theory or application, by where it enters | C's doubled, front halved |
| Sign of the operator on a cosine (+ for -) | application | C1, C2, C3, CG with the wrong sign; the front leads in the tough bands, f(0) > 0 |
| Slope signs reversed | theory if the coefficients are written with the wrong sign; application if a right formula is miscoded | S2, S3 with the wrong sign |
| Maximum of kappa3 taken from the sixteen table points | application | kappa3_max and s_kappa3_max off, tables right |
