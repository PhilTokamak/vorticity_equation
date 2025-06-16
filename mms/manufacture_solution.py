import sympy as sp

# Define symbols
x, y = sp.symbols('x y')

# Define ψ(x, y)
psi = (sp.sin(x))**3 * (sp.cos(5*y))**2

# Calculate ζ = -∇²ψ
psi_xx = sp.diff(psi, x, 2)
psi_yy = sp.diff(psi, y, 2)
zeta = (psi_xx + psi_yy)

# Calculate the first order derivate of ζ
zeta_x = sp.diff(zeta, x)
zeta_y = sp.diff(zeta, y)

# Calculate the first order derivate of ψ
psi_x = sp.diff(psi, x)
psi_y = sp.diff(psi, y)

# Calculate Jacobian: J(ζ, ψ) = ζ_x * ψ_y - ζ_y * ψ_x
jacobian = zeta_x * psi_y - zeta_y * psi_x

print("psi: ")
sp.pprint(psi)

print("\nzeta: ")
sp.pprint(zeta.simplify())

print("\nJacobian: ")
sp.pprint(jacobian.simplify())