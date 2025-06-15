// param.h
#ifndef PARAM_H
#define PARAM_H

#include <cmath>
#include <string>

// Space parameters
extern int N; // Use a number that is a exponential of 2 to avoid segmentation fault
              // This fault is probably due to multigrid solver for Poisson equation
extern double L;

extern double dx; // dx should not be read from param file;

// Time paramters
extern double dt;
extern int STEPS;
// Save data every 'SAVE_EVERY' steps
extern int SAVE_EVERY;

// Viscosity
extern double nu;

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

extern Scheme scheme;

// Read parameters from a param file
void load_parameters(const std::string& filename);


#endif // PARAM_H