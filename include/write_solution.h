#ifndef WRITE_SOLUTION_H
#define WRITE_SOLUTION_H

#include "grid.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void save_grid(const Grid& f, int step);
void save_vector(std::vector<double> vec, std::string filename);

#endif // WRITE_SOLUTION_H