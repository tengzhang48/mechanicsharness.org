# Separation between unequal strips — answer

## Results

| Quantity | Value |
|---|---|
| 1. Moment per strip $M$ | $M = 0.2126\ \mathrm{N\cdot m}$ ($0.21263$) |
| 2. Advance of front beyond starter-cut root | $L = 10.80\ \mathrm{mm}$ ($10.800$) |
| 3. Opening at mid-bridged span | $\delta(L/2) = 38.2\ \mu\mathrm{m}$ ($38.22$) |

The instant asked for is when the connector at the **original starter-cut root** first reaches
$\delta_c$ (opening grows monotonically from the front, so the root connector is the first to
become inactive). At that instant the bridged zone spans exactly $0<x<L$, with
$\delta(0)=0$, $\delta(L)=\delta_c$.

## Mechanical model and boundary conditions

Work per unit width with Euler–Bernoulli arms (axial stretch and shear negligible; no axial
force since there are no end forces):

- $d_1 = E_1h_1^3/12 = 1.9800\ \mathrm{N\cdot m}$, $d_2 = E_2h_2^3/12 = 11.5425\ \mathrm{N\cdot m}$,
  combined $\bar d = (1/d_1+1/d_2)^{-1} = 1.690083\ \mathrm{N\cdot m}$.
- The applied couples are self-equilibrating (zero net force and moment), so the far-ahead bonded
  bilayer stays straight and stress free: both arms are **clamped at the advancing front**,
  $\delta(0)=0$, $\delta'(0)=0$.
- On the bridged span $0<x<L$ the closing traction $t=\sigma_b(1-\delta/\delta_c)$ acts on both
  arms, giving (with $u\equiv 1-\delta/\delta_c$)
  $$u''''=\lambda^4 u,\qquad \lambda=\Big(\frac{\sigma_b}{\delta_c\bar d}\Big)^{1/4}=155.964\ \mathrm{m^{-1}}\quad(1/\lambda=6.412\ \mathrm{mm}).$$
- **Front condition**: within beam theory the energy release rate at the front is set by the tip
  moments $m_1(0)=-m_2(0)=\bar d\,\delta''(0)$:
  $$G_{tip}=\frac{m_1(0)^2}{2d_1}+\frac{m_2(0)^2}{2d_2}=\frac{(\bar d\,\delta''(0))^2}{2\bar d}=\Gamma_0
  \;\Rightarrow\; u''(0)=-\frac{1}{\delta_c}\sqrt{\frac{2\Gamma_0}{\bar d}}.$$
- **Trailing edge (root) conditions**: $\delta(L)=\delta_c$ (defines the instant), and zero shear
  at $x=L$ ($u'''(L)=0$): the traction vanishes at $\delta_c$ and the starter-cut region $L<x<a$
  carries no load, so shear is zero and moment constant ($=\pm m$) all the way to the free ends,
  where couples $\pm M$ act with no forces. (The starter-cut length 20 mm does not otherwise
  enter.)
- Quasi-static monotone loading: $\delta$ increases everywhere in time, so all connectors are still
  on the linear softening branch; no unloading occurs before the first one reaches $\delta_c$.

## Solution

$u=A\cosh\lambda x+B\sinh\lambda x+C\cos\lambda x-B\sin\lambda x$, with
$A=\tfrac{1-k}{2},\;C=\tfrac{1+k}{2},\;k=\sqrt{2\Gamma_0/(\sigma_b\delta_c)}=\sqrt{0.07}=0.264575$
(from the conditions at $x=0$). The two conditions at $x=L$ give a single equation for
$s=\lambda L$:
$$\cosh s\cos s+1=k\sinh s\sin s\;\Rightarrow\;s=1.684422 .$$

Hence $L=s/\lambda=10.800\ \mathrm{mm}$. The applied end moment per unit width is the moment at
$x=L$, $m=\bar d\,\delta''(L)$, which by the J-integral equals
$\sqrt{2\bar d(\Gamma_0+\tfrac12\sigma_b\delta_c)}=8.50503\ \mathrm{N}$, so
$$M=mb=8.50503\times0.025=\mathbf{0.21263\ N\cdot m}.$$
At mid-span $u(L/2)=0.80891$, i.e. $\delta(L/2)=\delta_c(1-u)=\mathbf{38.22\ \mu m}$.

## Checks

- **J-integral (energy) balance**: $m^2/2\bar d = 21.4000000\ \mathrm{J/m^2}$ vs
  $\Gamma_0+\int_0^{\delta_c}t\,d\delta = 1.40+0.5(0.2\mathrm{MPa})(0.2\mathrm{mm})=21.400\ \mathrm{J/m^2}$ — agree to $1.5\times10^{-15}$ relative.
- **Statics of an arm**: $m-m_{tip}=6.32966\ \mathrm{N}$ vs $\int_0^L t(x)\,x\,dx=6.32966\ \mathrm{N}$ (moment of bridging tractions about the front) — agree to $8\times10^{-15}$.
- The computed profile is monotone, $0\le\delta\le\delta_c$ on $[0,L]$, $t(0)=\sigma_b$, $t(L)=0$;
  slopes stay small ($\delta'\le0.043$, end rotation $\approx0.14$ rad), consistent with small-rotation
  beam theory.
- Sensible limits: $k\to0$ (no bridging, $\sigma_b\delta_c\to\infty$) gives $s\to$ root of
  $\cosh s\cos s=-1$, $L\to0$; with $\Gamma_0\to0$, $s\to\pi/2$.

Numbers to quote: **M = 0.2126 N·m per strip, L = 10.80 mm, δ(L/2) = 38.2 µm.**
