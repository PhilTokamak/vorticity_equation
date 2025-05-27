#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "grid.h"
#include "general_operatos.h"

void smooth(Grid& psi, const Grid& rhs, int iterations, double h);
Grid residual(const Grid& psi, const Grid& rhs, double h);
Grid restrict_grid(const Grid& fine);
Grid prolong_grid(const Grid& coarse);
Grid multigrid_vcycle(const Grid& rhs, int level, double h);
inline Grid multigrid_solve(const Grid& rhs, double h) {
    int levels = int(log2(rhs.size())) - 2;
    return multigrid_vcycle(rhs, levels, h);
}

#endif // MULTIGRID_H