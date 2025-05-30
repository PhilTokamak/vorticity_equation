# 2D Vorticity Equation
The repo contains the solver for vorticity equations in 2D, it supports several different schemes to discretize the Jacobian term.

## Build process
1. Clone this repo
2. Run the commands below
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

All the simulation parameters are in `include/param.hpp`. After changing any parameter, run `make` again in *build* directory to rebuild the project, and restart the simulation.

## OpenMP related
This project is parallelized with OpenMP, the default setting is to OpenMP. Please check the output aftering running `cmake ..` in *build* directory to see if OpenMP is enabled. In order to use OpenMP to accelerate the simulation, run, for instance,
```bash
export OMP_NUM_THREADS=10
```
to set the number of OpenMP threads to 10, before starting the simulation.

To disable OpenMP, please run
```bash
make clean
cmake -DUSE_OPENMP=OFF ..
```
in *build* directory.

#### MacOS
In MacOS, the default Apple clang compiler doesn't support OpenMP, in order to use it on MacOS, please use *homebrew* to install GNU gcc by running
```bash
brew install gcc
```
then run
```bash
which g++-14
```
to print the path for the GNU compiler, and either run the commands below or added them to `~/.zshrc` or `~/.bashrc` to enable the GNU compiler:
```bash
export CC=/opt/homebrew/bin/gcc-14
export CXX=/opt/homebrew/bin/g++-14
```
The exact path depends on the output from `which g++-14` command.