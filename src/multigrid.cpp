#include "multigrid.h"

/**
 * Using Gauss-Seidel method as pre- and post-smoothing to reduce high frequency errors.
 */
void smooth(Grid& psi, const Grid& rhs, int iterations, double h) {
    int n = psi.size();
    Grid psi_new = zero_grid(n);
    for (int it = 0; it < iterations; ++it) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                psi_new[i][j] = 0.25 * (psi[idx(i+1, n)][j] + psi[idx(i-1, n)][j]
                            + psi[i][idx(j+1, n)] + psi[i][idx(j-1, n)]
                            - h * h * rhs[i][j]);
            }
        }
        psi = psi_new;
    }
}

/**
 * Compute residual errors ：res = rhs - Lψ
 */
Grid residual(const Grid& psi, const Grid& rhs, double h) {
    return rhs - laplacian(psi, h);
}

/**
 * Downsampling (restricting) the residual error to a coarser grid.
 */
Grid restrict_grid(const Grid& fine) {
    int Nc = fine.size() / 2;
    int n = fine.size();
    Grid coarse = zero_grid(Nc);
    for (int i = 0; i < Nc; ++i) {
        for (int j = 0; j < Nc; ++j) {
            coarse[i][j] = 0.25 * (fine[idx(2*i, n)][idx(2*j, n)]
                                    + fine[idx(2*i+1, n)][idx(2*j, n)]
                                    + fine[idx(2*i, n)][idx(2*j+1, n)]
                                    + fine[idx(2*i+1, n)][idx(2*j+1, n)]);
        }
    }
    return coarse;
}
/**
 * Interpolating (prolonging) a correction computed on a coarser grid into a finer grid.
 */
Grid prolong_grid(const Grid& coarse) {
    int Nc = coarse.size();
    int Nf = Nc * 2;
    Grid fine = zero_grid(Nf);
    for (int i = 0; i < Nc; ++i) {
        for (int j = 0; j < Nc; ++j) {
            for (int di = 0; di < 2; ++di) {
                for (int dj = 0; dj < 2; ++dj) {
                    fine[idx(2*i + di, Nf)][idx(2*j + dj, Nf)] = coarse[i][j];
                }
            }
        }
    }
    return fine;
}

/**
 *  V-cycle Multigrid main function
    - level = 0 means the coarsest gird
    - Otherwise, recursively downsampling, computing errors and prolonging
 */
Grid multigrid_vcycle(const Grid& rhs, int level, double h) {
    int n = rhs.size();
    Grid psi = zero_grid(n);

    // Stop recursion at smallest grid size, otherwise continue recursion
    if (n <= 4 || level == 0) {
        smooth(psi, rhs, 50, h);
        return psi;
    }

    // Pre-Smoothing
    smooth(psi, rhs, 3, h);
    // Compute Residual Errors
    Grid r = residual(psi, rhs, h);
    // Restriction
    Grid r_coarse = restrict_grid(r);
    // Recursively compute errors
    Grid e_coarse = multigrid_vcycle(r_coarse, level - 1, 2 * h);
    // Prolongation
    Grid e_fine = prolong_grid(e_coarse);
    // Correction
    psi = psi + e_fine;
    // Post-Smoothing
    smooth(psi, rhs, 3, h);

    return psi;
}