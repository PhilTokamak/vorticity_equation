# 2D Vorticity Equation
The repo contains the solver for vorticity equations in 2D, it supports several different schemes to discretize the Jacobian term.

## Build process
1. Clone this repo
2. Run the commands below to use default compiler on your machine
```bash
cd vorticity_equation
mkdir build
cd build
cmake ..
make
```

### Build with non-default compilers
Currently supported compilers include GNU gcc compiler and Apple clang compiler. These compilers can be chosen using cmake presets. To do so, especially on macOS, after running `cd build` command, run
```bash
cmake --list-presets=configure ..
```
to check the available configure presets.

For the installation of GNU gcc on macOS, please refer to [this section](#macos). Depending on the version of GNU `gcc` on your machine, for instance, if you are using a version of 15.x, run the commands below to set the compiler to `g++-15`:
```bash
cmake --preset=g++-15 ..
cmake ..
```
Finally run
```bash
make
```
to build the project.

## Start a simulation
To start a simulation, just run
```bash
./start_simulation.sh
```

All the simulation parameters should be included in a `params_in.txt` file in build folder, otherwise default parameters will be used.

An example of `params_in.txt` file is given as follows:
```
# Simulation parameters
// This is an example of param file
N = 128
L = 6.283185307179586
dt = 0.01
STEPS = 500
SAVE_EVERY = 100
nu = 0
scheme = ARAKAWA_2
```

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
in *build* directory. Or alternatively, you can use directly `cmake --preset=clang ..` to use Apple clang for building, as described in [this section](#build-with-non-default-compilers).

#### MacOS
In MacOS, the default Apple clang compiler doesn't support OpenMP, in order to use it on MacOS, please use *homebrew* to install GNU gcc by running
```bash
brew install gcc
```
then depending on the version of gcc you installed, run either
```bash
which g++-14
```
or
```bash
which g++-15
```
to print the path for the GNU compiler.

#### Other machines
After having found the path the compiler, either run the commands below or added them to `~/.zshrc` or `~/.bashrc` and run `source ~/.zshrc` or `source ~/.bashrc` to enable the GNU compiler:
```bash
export CC=/opt/homebrew/bin/gcc-15
export CXX=/opt/homebrew/bin/g++-15
```
Note that The exact path depends on the the compiler you want to use.