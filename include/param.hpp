// param.hpp
#ifndef PARAM_HPP
#define PARAM_HPP

#include <cmath>

//Space parameters
constexpr int N = 128; // Use a number that is a exponential of 2 to avoid segmentation fault
                       // This fault is probably due to multigrid solver for Poisson equation
constexpr double L = 2 * M_PI;

// Time paramters
constexpr double dt = 0.05;
constexpr int STEPS = 4000;
// Save data every 'SAVE_EVERY' steps
constexpr int SAVE_EVERY = 100;

// Viscosity
constexpr double nu = 1e-4;

// Scheme choice
constexpr bool use_arakawa = true;


#endif // PARAM_HPP