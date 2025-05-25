#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

constexpr int N = 64;
constexpr double L = 2 * M_PI;
constexpr double dx = L / N;
constexpr double dt = 0.05;
constexpr double nu = 1.0e-4;
constexpr int STEPS = 10000;
constexpr int SAVE_EVERY = 100;

using Grid = std::vector<std::vector<double>>;

Grid zero_grid(int n) {
    return Grid(n, std::vector<double>(n, 0.0));
}

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

int idx(int i, int n) {
    return (i + n) % n;
}

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

Grid arakawa_jacobian(const Grid& zeta, const Grid& psi) {
    auto n = zeta.size();
    Grid J = zero_grid(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double jplusplus = (zeta[idx(i + 1, n)][j] - zeta[idx(i - 1, n)][j])
                                * (psi[i][idx(j + 1, n)] - psi[i][idx(j - 1, n)])
                                - (zeta[i][idx(j + 1, n)] - zeta[i][idx(j - 1, n)])
                                * (psi[idx(i + 1, n)][j] - psi[idx(i - 1, n)][j]);

            double jpluscross = zeta[idx(i + 1, n)][j]
                                * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                    - psi[idx(i + 1, n)][idx(j - 1, n)])
                                - zeta[idx(i - 1, n)][j]
                                * (psi[idx(i - 1, n)][idx(j + 1, n)]
                                    - psi[idx(i - 1, n)][idx(j - 1, n)])
                                - zeta[i][idx(j + 1, n)]
                                * (psi[idx(i + 1, n)][idx(j + 1, n)]
                                    - psi[idx(i - 1, n)][idx(j + 1, n)])
                                + zeta[i][idx(j - 1, n)]
                                * (psi[idx(i + 1, n)][idx(j - 1, n)]
                                    - psi[idx(i - 1, n)][idx(j - 1, n)]);

            double jcrossplus = zeta[idx(i + 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 1, n)] - psi[idx(i + 1, n)][j])
                                - zeta[idx(i - 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i - 1, n)][j] - psi[i][idx(j - 1, n)])
                                - zeta[idx(i - 1, n)][idx(j + 1, n)]
                                    * (psi[i][idx(j + 1, n)] - psi[idx(i - 1, n)][j])
                                + zeta[idx(i + 1, n)][idx(j - 1, n)]
                                    * (psi[idx(i + 1, n)][j] - psi[i][idx(j - 1, n)]);

            J[i][j] = (jplusplus + jpluscross + jcrossplus) / (12 * dx * dx);
        }
    }
    return J;
}

void smooth(Grid& psi, const Grid& rhs, int iterations, double h) {
    int n = psi.size();
    Grid psi_new = zero_grid(n);
    for (int it = 0; it < iterations; ++it) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                psi_new[i][j] = 0.25 * (psi[idx(i+1, n)][j] + psi[idx(i-1, n)][j]
                            + psi[i][idx(j+1, n)] + psi[i][idx(j-1, n)]
                            - h * h * rhs[i][j]);
            }
        }
        psi = psi_new;
    }
}

Grid residual(const Grid& psi, const Grid& rhs, double h) {
    return rhs - laplacian(psi, h);
}

Grid restrict_grid(const Grid& fine) {
    int Nc = fine.size() / 2;
    int n = fine.size();
    Grid coarse = zero_grid(Nc);
    for (int i = 0; i < Nc; ++i) {
        for (int j = 0; j < Nc; ++j) {
            coarse[i][j] = 0.25 * (fine[idx(2*i, n)][idx(2*j, n)] + fine[idx(2*i+1, n)][idx(2*j, n)]
                                    + fine[idx(2*i, n)][idx(2*j+1, n)] + fine[idx(2*i+1, n)][idx(2*j+1, n)]);
        }
    }
    return coarse;
}

Grid prolong_grid(const Grid& coarse) {
    int Nc = coarse.size();
    int Nf = Nc * 2;
    Grid fine = zero_grid(Nf);
    for (int i = 0; i < Nc; ++i) {
        for (int j = 0; j < Nc; ++j) {
            for (int di = 0; di < 2; ++di) {
                for (int dj = 0; dj < 2; ++dj) {
                    fine[idx(2*i + di, Nf)][idx(2*j + dj, Nf)] = coarse[i][j];
                }
            }
        }
    }
    return fine;
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

    psi = psi + e_fine;

    smooth(psi, rhs, 3, h);
    return psi;
}

Grid multigrid_solve(const Grid& rhs, double h) {
    int levels = int(log2(rhs.size())) - 2;
    return multigrid_vcycle(rhs, levels, h);
}


Grid rhs(const Grid& zeta) {
    Grid psi = multigrid_solve(zeta, dx);
    Grid jac = arakawa_jacobian(zeta, psi);
    Grid diff = nu * laplacian(zeta, dx);
    return jac + diff;
}

void rk4_step(Grid& zeta) {
    Grid k1 = rhs(zeta);
    Grid k2 = rhs(zeta + 0.5 * dt * k1);
    Grid k3 = rhs(zeta + 0.5 * dt * k2);
    Grid k4 = rhs(zeta + dt * k3);
    zeta = zeta + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void initialize(Grid& zeta) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * dx;
            double y = j * dx;
            zeta[i][j] = std::exp(-10.0 * ((x - 3 * M_PI / 4) * (x - 3 * M_PI / 4)
                                            + (y - M_PI) * (y - M_PI)))
                        + std::exp(-10.0 * ((x - 5 * M_PI / 4) * (x - 5 * M_PI / 4)
                                            + (y - M_PI) * (y - M_PI)));
        }
    }
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
    Grid zeta = zero_grid(N);
    initialize(zeta);
    system("rm -rf output");
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