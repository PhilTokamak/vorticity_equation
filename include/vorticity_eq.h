#ifndef VORTICITY_EQ_H
#define VORTICITY_EQ_H

#include "grid.hpp"
#include "multigrid.h"

constexpr double dt = 0.5;
constexpr double nu = 0;
constexpr int STEPS = 4000;
constexpr int SAVE_EVERY = 100;

Grid rhs(const Grid& zeta);
void rk4_step(Grid& zeta);
void initialize(Grid& zeta);

#endif // VORTICITY_EQ_H