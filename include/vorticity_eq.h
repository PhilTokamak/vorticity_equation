#ifndef VORTICITY_EQ_H
#define VORTICITY_EQ_H

#include "param.hpp"
#include "grid.h"
#include "multigrid.h"

Grid rhs(const Grid& zeta);
void rk4_step(Grid& zeta);
void initialize(Grid& zeta);

#endif // VORTICITY_EQ_H