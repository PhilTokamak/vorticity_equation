#ifndef VORTICITY_EQ_H
#define VORTICITY_EQ_H

#include "param.h"
#include "grid.h"
#include "multigrid.h"
#include "utility.h"
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

Grid fd_jacobian(const Grid& zeta, const Grid& psi);
Grid arakawa_jacobian(const Grid& zeta, const Grid& psi);
Grid arakawa_jacobian_4th_order(const Grid& zeta, const Grid& psi);
Grid rhs(const Grid& zeta);
void rk4_step(Grid& zeta);
void initialize(Grid& zeta);

#endif // VORTICITY_EQ_H