# A noninteger hardening exponent

A stationary, sharp crack occupies the negative x-axis in a homogeneous body in **plane strain**, with traction-free faces at θ = ±π. Polar coordinates are centered at its tip; θ = 0 points ahead of it. Loading is symmetric opening, increases proportionally, and has positive energy flux J. Use small-strain kinematics and the following idealized, fully plastic total-deformation law:

\[
\epsilon_{ij}=\frac32\alpha\epsilon_0
 \left(\frac{\sigma_e}{\sigma_0}\right)^{n-1}\frac{s_{ij}}{\sigma_0},
\qquad s_{ij}=\sigma_{ij}-\frac13\sigma_{kk}\delta_{ij},
\qquad \sigma_e=\sqrt{\frac32s_{ij}s_{ij}}.
\]

Here n = **6.5**, α = 0.02, σ₀ = 400 MPa, E = 200 GPa, and ε₀ = σ₀/E. E only sets the reference strain; do not add an elastic strain term to this limiting law. Indices run over all three spatial directions, repeated indices are summed, and εᵣθ is tensor shear strain. Plane strain means εzz = εrz = εθz = 0.

Let W = n σᵢⱼεᵢⱼ/(n+1) be strain energy per volume. On a counterclockwise circular contour, with outward radial normal N and traction tᵢ = σᵢⱼNⱼ, define

\[
J=\int_{-\pi}^{\pi}\left(W\cos\theta-t_i\frac{\partial u_i}{\partial x}\right)r\,d\theta.
\]

Crack-face contributions vanish. J has units J/m². Let U(R) = ∫₀ᴿ∫₋π^π W r dθ dr, the energy inside radius R **per unit out-of-plane thickness**, in J/m. Work with the leading tip field and take R within its idealized range.

Compute these three dimensionless numbers:

1. p in σθθ(r,0) ∝ r⁻ᵖ, to absolute tolerance **10⁻⁵**.
2. H defined by σθθ(r,0) = H σ₀ [J/(α ε₀ σ₀ r)]ᵖ, to absolute tolerance **0.005**.
3. C = U(R)/(J R), to absolute tolerance **0.002**.

Submit the numbers, a brief method, and one verification. No particular code format is required; allow about 45 minutes using NumPy/SciPy.

Computational identities, if useful: an Airy function Φ gives σrr = Φ,r/r + Φ,θθ/r², σθθ = Φ,rr, σrθ = −∂r(Φ,θ/r). In polar components, εrr = ur,r, εθθ = (ur + uθ,θ)/r, and 2εrθ = ur,θ/r + uθ,r − uθ/r. Cartesian differentiation of displacement must include rotation of the polar basis.
