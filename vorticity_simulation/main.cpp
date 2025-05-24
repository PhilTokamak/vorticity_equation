#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

constexpr int N = 128;
constexpr double L = 2 * M_PI;
constexpr double dx = L / N;
constexpr double dt = 0.01;
constexpr int STEPS = 2000;
constexpr int SAVE_EVERY = 10;

using Grid = std::vector<std::vector<double>>;

int idx(int i, int n = N) {
    return (i + n) % n;
}

Grid zero_grid(int n) {
    return Grid(n, std::vector<double>(n, 0.0));
}

void initialize(Grid& zeta) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double x = i * dx;
            double y = j * dx;
            zeta[i][j] = std::exp(-10.0 * ((x - M_PI) * (x - M_PI) + (y - M_PI) * (y - M_PI)));
        }
}

Grid restrict_grid(const Grid& fine) {
    int Nc = fine.size() / 2;
    Grid coarse = zero_grid(Nc);
    for (int i = 1; i < Nc - 1; ++i)
        for (int j = 1; j < Nc - 1; ++j)
            coarse[i][j] = 0.25 * (fine[2*i][2*j] + fine[2*i+1][2*j] + fine[2*i][2*j+1] + fine[2*i+1][2*j+1]);
    return coarse;
}

Grid prolong_grid(const Grid& coarse) {
    int Nc = coarse.size();
    int Nf = Nc * 2;
    Grid fine = zero_grid(Nf);
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j) {
            for (int di = 0; di < 2; ++di)
                for (int dj = 0; dj < 2; ++dj)
                    fine[2*i + di][2*j + dj] = coarse[i][j];
        }
    return fine;
}

Grid residual(const Grid& psi, const Grid& rhs, double h) {
    int n = psi.size();
    Grid r = zero_grid(n);
    for (int i = 1; i < n - 1; ++i)
        for (int j = 1; j < n - 1; ++j)
            r[i][j] = rhs[i][j] - (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - 4 * psi[i][j]) / (h * h);
    return r;
}

void smooth(Grid& psi, const Grid& rhs, int iterations, double h) {
    int n = psi.size();
    for (int it = 0; it < iterations; ++it) {
        for (int i = 1; i < n - 1; ++i)
            for (int j = 1; j < n - 1; ++j)
                psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + h * h * rhs[i][j]);
    }
}

Grid multigrid_vcycle(const Grid& rhs, int level, double h) {
    int n = rhs.size();
    Grid psi = zero_grid(n);

    if (n <= 4 || level == 0) {
        smooth(psi, rhs, 50, h);
        return psi;
    }

    smooth(psi, rhs, 3, h);
    Grid r = residual(psi, rhs, h);
    Grid r_coarse = restrict_grid(r);
    Grid e_coarse = multigrid_vcycle(r_coarse, level - 1, 2 * h);
    Grid e_fine = prolong_grid(e_coarse);

    for (int i = 1; i < n - 1; ++i)
        for (int j = 1; j < n - 1; ++j)
            psi[i][j] += e_fine[i][j];

    smooth(psi, rhs, 3, h);
    return psi;
}

Grid multigrid_solve(const Grid& rhs, int levels) {
    return multigrid_vcycle(rhs, levels, 1.0 / rhs.size());
}

Grid arakawa_jacobian(const Grid& zeta, const Grid& psi) {
    Grid J(N, std::vector<double>(N, 0.0));
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            double jpp = (zeta[idx(i + 1)][j] - zeta[idx(i - 1)][j]) * (psi[i][idx(j + 1)] - psi[i][idx(j - 1)])
                       - (zeta[i][idx(j + 1)] - zeta[i][idx(j - 1)]) * (psi[idx(i + 1)][j] - psi[idx(i - 1)][j]);

            double jxp = zeta[idx(i + 1)][j] * (psi[idx(i + 1)][idx(j + 1)] - psi[idx(i + 1)][idx(j - 1)])
                       - zeta[idx(i - 1)][j] * (psi[idx(i - 1)][idx(j + 1)] - psi[idx(i - 1)][idx(j - 1)])
                       - zeta[i][idx(j + 1)] * (psi[idx(i + 1)][idx(j + 1)] - psi[idx(i - 1)][idx(j + 1)])
                       + zeta[i][idx(j - 1)] * (psi[idx(i + 1)][idx(j - 1)] - psi[idx(i - 1)][idx(j - 1)]);

            double jpx = zeta[idx(i + 1)][idx(j + 1)] * (psi[i][idx(j + 1)] - psi[idx(i + 1)][j])
                       - zeta[idx(i - 1)][idx(j + 1)] * (psi[i][idx(j + 1)] - psi[idx(i - 1)][j])
                       - zeta[idx(i + 1)][idx(j - 1)] * (psi[i][idx(j - 1)] - psi[idx(i + 1)][j])
                       + zeta[idx(i - 1)][idx(j - 1)] * (psi[i][idx(j - 1)] - psi[idx(i - 1)][j]);

            J[i][j] = (jpp + jxp + jpx) / (12 * dx * dx);
        }
    }
    return J;
}

void rk4_step(Grid& zeta) {
    auto k1 = arakawa_jacobian(zeta, multigrid_solve(zeta, 4));
    Grid tmp1 = zeta, tmp2 = zeta, tmp3 = zeta;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp1[i][j] -= 0.5 * dt * k1[i][j];

    auto k2 = arakawa_jacobian(tmp1, multigrid_solve(tmp1, 4));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp2[i][j] -= 0.5 * dt * k2[i][j];

    auto k3 = arakawa_jacobian(tmp2, multigrid_solve(tmp2, 4));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) tmp3[i][j] -= dt * k3[i][j];

    auto k4 = arakawa_jacobian(tmp3, multigrid_solve(tmp3, 4));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            zeta[i][j] += dt * (k1[i][j] + 2 * k2[i][j] + 2 * k3[i][j] + k4[i][j]) / 6.0;
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
            std::cout << "Step " << step << std::endl;
        }
        rk4_step(zeta);
    }
    std::cout << "Simulation complete.\n";
    return 0;
}