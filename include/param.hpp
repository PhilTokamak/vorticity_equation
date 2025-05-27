// param.hpp
#ifndef PARAM_HPP
#define PARAM_HPP

#include <cmath>

//Space parameters
constexpr int N = 256; // Use a number that is a exponential of 2 to avoid segmentation fault
                       // This fault is probably due to multigrid solver for Poisson equation
constexpr double L = 2 * M_PI;

// Time paramters
constexpr double dt = 0.05;
constexpr int STEPS = 4000;
// Save data every 'SAVE_EVERY' steps
constexpr int SAVE_EVERY = 100;

// Viscosity
constexpr double nu = 0;

/**
 *  Scheme choice for Jacobian:
 *      1. ARAKAWA_2: 2nd order nine-point Arakawa scheme
 *      2. ARAKAWA_4: 4th order thirteen-point Arakawa scheme
 *      3. CENTERED_2: 2nd order centered scheme
*/
enum Scheme {
    ARAKAWA_2 = 1,
    ARAKAWA_4 = 2,
    CENTERED_2 = 3,
};

constexpr Scheme scheme = ARAKAWA_2;


#endif // PARAM_HPP