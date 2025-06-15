#ifndef UTILITY_H
#define UTILITY_H

#include "grid.h"
#include <random>

std::string trim_spaces(const std::string& line);
Grid laplacian(const Grid& f, double h);
double rand_real(double lower_bound = 0.0, double upper_bound = 1.0);

#endif // UTILITY_H