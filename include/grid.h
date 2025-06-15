#ifndef GRID_H
#define GRID_H

#include "param.h"
#include <vector>

class Grid {
private:
    int _n;
    std::vector<std::vector<double>> _data;

public:
    // Constructor, create n x n zero matrix
    Grid(int n) : _n(n), _data(n, std::vector<double>(n, 0.0)) {}
    int size() const { return _n; }

    // Overload [] operator to access the data
    std::vector<double>& operator[](int i) { return _data[i]; }
    const std::vector<double>& operator[](int i) const { return _data[i]; }

    // Overload operators for Grid

    // Addition between two Grids
    friend Grid operator+(const Grid& a, const Grid& b);
    // Subtraction between two Grids
    friend Grid operator-(const Grid& a, const Grid& b);
    // Multiply a Grid by a real number
    friend Grid operator*(double scalar, const Grid& g);
    // Element-wise multiplication between tow Grids
    friend Grid operator*(const Grid& a, const Grid& g);

    //Operations on Grid

    // Compute the sum of all elements in the grid
    double sum() const;
};

// Create a all-zero grid
inline Grid zero_grid(int n) {
    return Grid(n);
}

// Return periodic index
inline int idx(int i, int n) {
    return (i + n) % n;
}

#endif // GRID_H