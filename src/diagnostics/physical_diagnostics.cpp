#include "param.hpp"
#include "physical_diagnostics.h"

double compute_mean(const Grid& f) {
    return f.sum() * dx * dx / (L * L);
}

std::pair<Grid,Grid> compute_velocity(const Grid& psi, double h) {
    int n = psi.size();
    Grid vx(n), vy(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            vx[i][j] =   (psi[i][idx(j + 1, n)] - psi[i][idx(j - 1, n)]) / (2 * h);
            vy[i][j] = - (psi[idx(i + 1, n)][j] - psi[idx(i - 1, n)][j]) / (2 * h);
        }
    }
    return std::make_pair(vx, vy);
}