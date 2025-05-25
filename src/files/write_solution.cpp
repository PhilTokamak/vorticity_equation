#include "write_solution.h"

void save_csv(const Grid& f, int step) {
    std::ostringstream filename;
    filename << "output/zeta_" << step << ".csv";
    std::ofstream file(filename.str());
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            file << f[i][j];
            if (i < N - 1) file << ",";
        }
        file << "\n";
    }
}