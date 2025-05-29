#include "vorticity_eq.h"
#include <stdexcept>

/**
 * Jacobian J(zeta, psi) discretized by second order finite difference
*/
Grid fd_jacobian(const Grid& zeta, const Grid& psi) {
    auto n = zeta.size();
    Grid J = zero_grid(n);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(zeta, psi, J, n)
#endif //_OPENMP
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            J[i][j] = ((zeta[idx(i+1, n)][j] - zeta[idx(i-1, n)][j])
                        * (psi[i][idx(j+1, n)] - psi[i][idx(j-1, n)])
                        - (zeta[i][idx(j+1, n)] - zeta[i][idx(j-1, n)])
                        * (psi[idx(i+1, n)][j] - psi[idx(i-1, n)][j]))
                      / (4 * dx * dx);
        }
    }
    return J;
}

/**
 * Arakawa Jacobian (2nd order nine-point scheme)
 * Eq. (36) ~ (40), (44), Arakawa_1966
 */
Grid arakawa_jacobian(const Grid& zeta, const Grid& psi) {
    auto n = zeta.size();
    Grid J = zero_grid(n);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(zeta, psi, J, n)
#endif //_OPENMP
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double jplusplus = (zeta[idx(i + 1, n)][j] - zeta[idx(i - 1, n)][j])
                                * (psi[i][idx(j + 1, n)] - psi[i][idx(j - 1, n)])
                                - (zeta[i][idx(j + 1, n)] - zeta[i][idx(j - 1, n)])
                                * (psi[idx(i + 1, n)][j] - psi[idx(i - 1, n)][j]);

            double jpluscross = zeta[idx(i + 1, n)][j]
                                    * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                        - psi[idx(i + 1, n)][idx(j - 1, n)])
                                - zeta[idx(i - 1, n)][j]
                                    * (psi[idx(i - 1, n)][idx(j + 1, n)]
                                        - psi[idx(i - 1, n)][idx(j - 1, n)])
                                - zeta[i][idx(j + 1, n)]
                                    * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                        - psi[idx(i - 1, n)][idx(j + 1, n)])
                                + zeta[i][idx(j - 1, n)]
                                    * (psi[idx(i + 1, n)][idx(j - 1, n)]
                                        - psi[idx(i - 1, n)][idx(j - 1, n)]);

            double jcrossplus = zeta[idx(i + 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 1, n)] - psi[idx(i + 1, n)][j])
                                - zeta[idx(i - 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i - 1, n)][j] - psi[i][idx(j - 1, n)])
                                - zeta[idx(i - 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 1, n)] - psi[idx(i - 1, n)][j])
                                + zeta[idx(i + 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i + 1, n)][j] - psi[i][idx(j - 1, n)]);

            J[i][j] = (jplusplus + jpluscross + jcrossplus) / (12 * dx * dx);
        }
    }
    return J;
}

/**
 * Arakawa Jacobian (4th order thirteen-point scheme)
 * Eq. (39), (56) ~ (58), Arakawa_1966
 */
Grid arakawa_jacobian_4th_order(const Grid& zeta, const Grid& psi) {
    auto n = zeta.size();
    Grid J = zero_grid(n);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(zeta, psi, J, n)
#endif //_OPENMP
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double jcrosscross = (zeta[idx(i + 1, n)][idx(j + 1, n)]
                                    - zeta[idx(i - 1, n)][idx(j - 1, n)])
                                * (psi[idx(i - 1, n)][idx(j + 1, n)]
                                    - psi[idx(i + 1, n)][idx(j - 1, n)])
                                - (zeta[idx(i - 1, n)][idx(j + 1, n)]
                                    - zeta[idx(i + 1, n)][idx(j - 1, n)])
                                * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                    - psi[idx(i - 1, n)][idx(j - 1, n)]);

            double jcrossplus = zeta[idx(i + 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 2, n)] - psi[idx(i + 2, n)][j])
                                - zeta[idx(i - 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i - 2, n)][j] - psi[i][idx(j - 2, n)])
                                - zeta[idx(i - 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 2, n)] - psi[idx(i - 2, n)][j])
                                + zeta[idx(i + 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i + 2, n)][j] - psi[i][idx(j - 2, n)]);

            double jpluscross = zeta[idx(i + 2, n)][j]
                                    * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                        - psi[idx(i + 1, n)][idx(j - 1, n)])
                                - zeta[idx(i - 2, n)][j]
                                    * (psi[idx(i - 1, n)][idx(j + 1, n)]
                                        - psi[idx(i - 1, n)][idx(j - 1, n)])
                                - zeta[i][idx(j + 2, n)]
                                    * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                        - psi[idx(i - 1, n)][idx(j + 1, n)])
                                + zeta[i][idx(j - 2, n)]
                                    * (psi[idx(i + 1, n)][idx(j - 1, n)]
                                        - psi[idx(i - 1, n)][idx(j - 1, n)]);

            J[i][j] = (jcrosscross + jcrossplus + jpluscross) / (24 * dx * dx);
        }
    }
    return J;
}

Grid rhs(const Grid& zeta) {
    Grid psi = multigrid_solve(zeta, dx);

    Grid jac(psi.size());
    switch (scheme) {
        case ARAKAWA_2:
            jac = arakawa_jacobian(zeta, psi);
            break;
        case ARAKAWA_4:
            jac = arakawa_jacobian_4th_order(zeta, psi);
            break;
        case CENTERED_2:
            jac = fd_jacobian(zeta, psi);
            break;
        default:
            throw std::invalid_argument("Invalid scheme: "
                                    + std::to_string(static_cast<int>(scheme)));
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
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(zeta)
#endif //_OPENMP
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * dx;
            double y = j * dx;
            // Two Gaussian vortices
            zeta[i][j] = std::exp(-10.0 * ((x - 3.0 / 4.0 * L / 2.0) * (x - 3.0 / 4.0 * L / 2.0)
                                            + (y - L / 2.0) * (y - L / 2.0)))
                        + std::exp(-10.0 * ((x - 5.0 / 4.0 * L / 2.0) * (x - 5.0 / 4.0 * L / 2.0)
                                            + (y - L / 2.0) * (y - L / 2.0)));
            // Taylor-Green vortex
            // zeta[i][j] = std::cos(2 * M_PI * x / L) * std::cos(2 * M_PI * y / L);
        }
    }
}