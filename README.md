# 2D Vorticity Equation
The repo contains the solver for vorticity equations in 2D, it supports several different schemes to discretize the Jacobian term.

## Build process
1. Clone this repo
2. run the commands below
```bash
cd vorticity_equation
mkdir build
cd build
cmake ..
make
```

## Start a simulation
To start a simulation, just run
```bash
./start_simulation.sh
```

All the simulation parameters are in `include/param.hpp`. After changing any parameter, run `make` again in `build` directory to rebuild the project, and restart the simulation.