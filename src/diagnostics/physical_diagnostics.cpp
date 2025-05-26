#include "param.hpp"
#include "physical_diagnostics.h"

double compute_mean(const Grid& f) {
    return f.sum() * dx * dx / (L * L);
}