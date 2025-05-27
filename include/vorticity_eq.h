#ifndef VORTICITY_EQ_H
#define VORTICITY_EQ_H

#include "param.hpp"
#include "grid.h"
#include "multigrid.h"
#include "general_operatos.h"


Grid fd_jacobian(const Grid& zeta, const Grid& psi);
Grid arakawa_jacobian(const Grid& zeta, const Grid& psi);
Grid arakawa_jacobian_4th_order(const Grid& zeta, const Grid& psi);
Grid rhs(const Grid& zeta);
void rk4_step(Grid& zeta);
void initialize(Grid& zeta);

#endif // VORTICITY_EQ_H