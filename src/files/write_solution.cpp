#include "write_solution.h"

void save_grid(const Grid& f, int step, std::string grid_name) {
    std::ostringstream filename;
    filename << "output/" << grid_name << "_" << step << ".csv";
    std::ofstream file(filename.str());
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            file << f[i][j];
            if (i < N - 1) file << ",";
        }
        file << "\n";
    }
}

void save_vector(const std::vector<double>& vec, const std::string filename) {
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

std::string scheme_to_string(Scheme scheme) {
    switch (scheme) {
        case ARAKAWA_2: return "2nd order Arakawa";
        case ARAKAWA_4: return "4th order Arakawa";
        case CENTERED_2: return "2nd order Centered Difference";
        default: return "UNKNOWN";
    }
}

void print_parameters_to_screen() {
    std::cout << "Parameters:" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "L = " << L << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "Total steps = " << STEPS << std::endl;
    std::cout << "SAVE_EVERY = " << SAVE_EVERY << std::endl;
    std::cout << "Viscosity nu = " << nu << std::endl;
    std::cout << "Use " << scheme_to_string(scheme)
                << " Jacobian" << std::endl <<std::endl;
}

void save_parameters(const std::string filename) {
    std::ofstream out(filename);

    if (!out) {
        std::cerr << "Failed to open" << filename << " for writing!" << std::endl;
        return;
    }
    out << "Parameters:" << std::endl;
    out << "N = " << N << std::endl;
    out << "L = " << L << std::endl;
    out << "dt = " << dt << std::endl;
    out << "Total steps = " << STEPS << std::endl;
    out << "SAVE_EVERY = " << SAVE_EVERY << std::endl;
    out << "Viscosity nu = " << nu << std::endl;
    out << "Use " << scheme_to_string(scheme)
            << "Jacobian" << std::endl;
    out.close();
}