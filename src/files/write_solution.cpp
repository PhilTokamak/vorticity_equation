#include "write_solution.h"

void save_grid(const Grid& f, int step) {
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

void save_vector(std::vector<double> vec, std::string filename) {
    std::ofstream out(filename);

    if (!out) {
        std::cerr << "Failed to open" << filename << " for writing!" << std::endl;
        return;
    }

    int n = vec.size();
    for (int i = 0; i < n; ++i) {
        out << vec[i];
            if (i < n - 1) {
                out << ",";
            } else {
                out << std::endl;
            }
    }
    out.close();
}