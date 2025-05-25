#include "grid.hpp"
#include "multigrid.h"
#include "vorticity_eq.h"
#include "write_solution.h"

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