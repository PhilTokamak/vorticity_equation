#include "grid.h"
#include "multigrid.h"
#include "vorticity_eq.h"
#include "write_solution.h"
#include "physical_diagnostics.h"
#include <fstream>
#include <iostream>

int main() {

#ifdef _OPENMP
    std::cout << "OpenMP is enabled. Running in multiple threads." << std::endl;
#pragma omp parallel
{
    #pragma omp master
    std::cout << "Number of threads in parallel region: " << omp_get_num_threads() << std::endl;
}
    std::cout << "Max threads: " << omp_get_max_threads() <<std::endl;
#else
    std::cout << "Warning: OpenMP not enabled. Running serially." << std::endl;
#endif //_OPENMP
    Grid zeta = zero_grid(N);
    std::vector<double> T{};
    std::vector<double> vorticity_mean{}, kinetic_E_mean{}, enstrophy_mean{};

    // Print info about this simulation
    std::cout << "Start simulation..." << std::endl << std::endl;
    print_parameters_to_screen();

    // create a new folder to store results
    system("rm -rf output");
    system("mkdir -p output");
    system("mkdir -p output/diagnostics");

    // Define file names
    std::string param_out_filename ="output/parameters_out.txt";
    std::string time_filename ="output/T.csv";
    std::string vorticity_mean_filename = "output/diagnostics/mean_vorticity.csv";
    std::string kinetic_E_mean_filename = "output/diagnostics/mean_kinetic_E.csv";
    std::string enstrophy_mean_filename = "output/diagnostics/mean_enstrophy.csv";

    initialize(zeta);

    for (int step = 0; step <= STEPS; ++step) {
        if (step % SAVE_EVERY == 0) {
            save_grid(zeta, step);

            // Solve Poisson equation to get stream function
            Grid psi = multigrid_solve(zeta, dx);

            // Compute velocity field
            auto [ux, uv] = compute_velocity(psi, dx);

            // Compute mean values of some physical quantities:
            // Mean vorticity
            double vorticity = compute_mean(zeta);
            vorticity_mean.push_back(vorticity);
            // Mean kinetic energy
            double kinetic_E = compute_mean(ux * ux + uv * uv);
            kinetic_E_mean.push_back(kinetic_E);
            // Mean enstrophy
            double enstrophy = compute_mean(zeta * zeta);
            enstrophy_mean.push_back(enstrophy);

            T.push_back(step * dt);
            std::cout << "Step " << step << std::endl;
        }
        rk4_step(zeta);
    }

    save_parameters(param_out_filename);
    save_vector(T, time_filename);
    save_vector(vorticity_mean, vorticity_mean_filename);
    save_vector(kinetic_E_mean, kinetic_E_mean_filename);
    save_vector(enstrophy_mean, enstrophy_mean_filename);
    std::cout << "Simulation complete.\n";
    return 0;
}