#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <cmath>

constexpr int N = 64;
constexpr double L = 2 * M_PI;
constexpr double dx = L / N;

using Grid = std::vector<std::vector<double>>;

// Create a all-zero grid
inline Grid zero_grid(int n) {
    return Grid(n, std::vector<double>(n, 0.0));
}

// Return periodic index
inline int idx(int i, int n) {
    return (i + n) % n;
}

// Reload operators for Grid
inline Grid operator+(const Grid& a, const Grid& b) {
    int n = a.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

inline Grid operator*(double scalar, const Grid& g) {
    int n = g.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = scalar * g[i][j];
        }
    }
    return result;
}

inline Grid operator-(const Grid& a, const Grid& b) {
    int n = a.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}

#endif // GRID_HPP