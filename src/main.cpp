#include "grid.h"
#include "multigrid.h"
#include "vorticity_eq.h"
#include "write_solution.h"
#include "physical_diagnostics.h"
#include <fstream>
#include <iostream>

int main() {
    Grid zeta = zero_grid(N);
    initialize(zeta);
    system("rm -rf output");
    system("mkdir -p output");
    system("mkdir -p output/diagnostics");

    // open diagnostics.csv file to write mean values data
    std::string diag_filename = "output/diagnostics/mean_vorticity.csv";
    std::ofstream diag_file(diag_filename);

    if (!diag_file) {
        std::cerr << "Failed to open" << diag_filename << " for writing!" << std::endl;
        return 1;
    }

    for (int step = 0; step <= STEPS; ++step) {
        if (step % SAVE_EVERY == 0) {
            save_csv(zeta, step);

            // Compute mean values of some physical quantities
            double zeta_mean = compute_mean(zeta);
            diag_file << zeta_mean;
            if (STEPS - step >= SAVE_EVERY) {
                diag_file << ",";
            } else {
                diag_file << std::endl;
            }

            std::cout << "Step " << step << std::endl;
        }
        rk4_step(zeta);
    }
    diag_file.close();
    std::cout << "Simulation complete.\n";
    return 0;
}