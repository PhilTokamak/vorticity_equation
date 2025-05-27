#ifndef WRITE_SOLUTION_H
#define WRITE_SOLUTION_H

#include "param.hpp"
#include "grid.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void save_grid(const Grid& f, int step);
void save_vector(const std::vector<double>& vec, const std::string filename);
std::string scheme_to_string(Scheme scheme);
void print_parameters_to_screen();
void save_parameters(const std::string filename);

#endif // WRITE_SOLUTION_H