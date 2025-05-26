#include "vorticity_eq.h"

Grid rhs(const Grid& zeta) {
    Grid psi = multigrid_solve(zeta, dx);

    Grid jac(psi.size());
    if (use_arakawa) {
        jac = arakawa_jacobian(zeta, psi);
    } else {
        jac = fd_jacobian(zeta, psi);
    }

    Grid diff = nu * laplacian(zeta, dx);
    return jac + diff;
}

void rk4_step(Grid& zeta) {
    Grid k1 = rhs(zeta);
    Grid k2 = rhs(zeta + 0.5 * dt * k1);
    Grid k3 = rhs(zeta + 0.5 * dt * k2);
    Grid k4 = rhs(zeta + dt * k3);
    zeta = zeta + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void initialize(Grid& zeta) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * dx;
            double y = j * dx;
            zeta[i][j] = std::exp(-10.0 * ((x - 3.0 / 4.0 * L / 2.0) * (x - 3.0 / 4.0 * L / 2.0)
                                            + (y - L / 2.0) * (y - L / 2.0)))
                        + std::exp(-10.0 * ((x - 5.0 / 4.0 * L / 2.0) * (x - 5.0 / 4.0 * L / 2.0)
                                            + (y - L / 2.0) * (y - L / 2.0)));
        }
    }
}