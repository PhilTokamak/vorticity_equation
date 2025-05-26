#include "grid.h"

Grid operator+(const Grid& a, const Grid& b) {
    int n = a.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

Grid operator-(const Grid& a, const Grid& b) {
    int n = a.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}

Grid operator*(double scalar, const Grid& g) {
    int n = g.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = scalar * g[i][j];
        }
    }
    return result;
}

Grid operator*(const Grid& a, const Grid& b) {
    int n = a.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = a[i][j] * b[i][j];
        }
    }
    return result;
}



double Grid::sum() const {
    int n = size();
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result += _data[i][j];
        }
    }
    return result;
}