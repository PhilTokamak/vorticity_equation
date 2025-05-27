#include "general_operatos.h"

/**
 * Define laplacian, use idx to preserve periodicity
*/
Grid laplacian(const Grid& f, double h) {
    int n = f.size();
    Grid result = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = (f[idx(i + 1, n)][j] + f[idx(i - 1, n)][j]
                            + f[i][idx(j + 1, n)] + f[i][idx(j - 1, n)]
                            - 4.0 * f[i][j]) / (h * h);
        }
    }
    return result;
}