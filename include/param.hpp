// param.hpp
#ifndef PARAM_HPP
#define PARAM_HPP

#include <cmath>

//Space parameters
constexpr int N = 64;
constexpr double L = 2 * M_PI;

// Time paramters
constexpr double dt = 0.05;
constexpr int STEPS = 4000;
// Save data
constexpr int SAVE_EVERY = 100;

// viscosity
constexpr double nu = 1e-4;


#endif // PARAM_HPP