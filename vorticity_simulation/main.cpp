
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

constexpr int N = 128;
constexpr double L = 2*M_PI;
constexpr double dx = L / N;
constexpr double dt = 0.01;
constexpr int STEPS = 1000;
constexpr int SAVE_EVERY = 10;

using Grid = std::vector<std::vector<double>>;

void initialize(Grid& zeta) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * dx - 0.5;
            double y = j * dx - 0.5;
            zeta[i][j] = std::exp(-10.0 * ((x - M_PI) * (x - M_PI) + (y - M_PI) * (y - M_PI)));
        }
    }
}

int idx(int i) {
    return (i + N) % N;
}

Grid laplacian(const Grid& f) {
    Grid result(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i][j] = (f[idx(i + 1)][j] + f[idx(i - 1)][j] +
                            f[i][idx(j + 1)] + f[i][idx(j - 1)] - 4 * f[i][j]) / (dx * dx);
    return result;
}

Grid poisson_solver(const Grid& zeta, int iterations = 1000) {
    Grid psi(N, std::vector<double>(N, 0.0));
    for (int it = 0; it < iterations; ++it) {
        Grid new_psi = psi;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                new_psi[i][j] = 0.25 * (psi[idx(i + 1)][j] + psi[idx(i - 1)][j] +
                                        psi[i][idx(j + 1)] + psi[i][idx(j - 1)] +
                                        dx * dx * zeta[i][j]);
        psi = new_psi;
    }
    return psi;
}

Grid arakawa_jacobian(const Grid& zeta, const Grid& psi) {
    Grid J(N, std::vector<double>(N, 0.0));
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            double jp = (zeta[idx(i + 1)][j] - zeta[idx(i - 1)][j]) * (psi[i][idx(j + 1)] - psi[i][idx(j - 1)])
                      - (zeta[i][idx(j + 1)] - zeta[i][idx(j - 1)]) * (psi[idx(i + 1)][j] - psi[idx(i - 1)][j]);
            J[i][j] = jp / (4 * dx * dx);
        }
    }
    return J;
}

void rk4_step(Grid& zeta) {
    auto k1 = arakawa_jacobian(zeta, poisson_solver(zeta));
    Grid tmp1 = zeta, tmp2 = zeta, tmp3 = zeta;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp1[i][j] -= 0.5 * dt * k1[i][j];

    auto k2 = arakawa_jacobian(tmp1, poisson_solver(tmp1));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp2[i][j] -= 0.5 * dt * k2[i][j];

    auto k3 = arakawa_jacobian(tmp2, poisson_solver(tmp2));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp3[i][j] -= dt * k3[i][j];

    auto k4 = arakawa_jacobian(tmp3, poisson_solver(tmp3));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            zeta[i][j] -= dt * (k1[i][j] + 2 * k2[i][j] + 2 * k3[i][j] + k4[i][j]) / 6.0;
}

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

int main() {
    Grid zeta(N, std::vector<double>(N));
    initialize(zeta);
    system("mkdir -p output");
    for (int step = 0; step <= STEPS; ++step) {

        if (step % SAVE_EVERY == 0) {
            save_csv(zeta, step);
            std::cout << "Step " << step <<std::endl;
        }
        rk4_step(zeta);
    }
    std::cout << "Simulation complete.\n";
    return 0;
}
