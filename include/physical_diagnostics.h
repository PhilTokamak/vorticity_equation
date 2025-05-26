#ifndef PHYSICAL_DIAGNOSTICS_H
#define PHYSICAL_DIAGNOSTICS_H

#include <grid.h>

double compute_mean(const Grid& f);
std::pair<Grid,Grid> compute_velocity(const Grid& psi, double h);

#endif // PHYSICAL_DIAGNOSTICS_H