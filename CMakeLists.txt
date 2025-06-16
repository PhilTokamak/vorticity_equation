
cmake_minimum_required(VERSION 3.10)

# Define project name
project(VorticitySimulation LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Adding option: whether to use OpenMP or not, default is to use OpenMP (ON);
# When using OpenMP, please run `export OMP_NUM_THREADS=10` to specfify the
# number of OpenMP threads;
# If not using OpenMP, use `cmake -DUSE_OPENMP=OFF ..` in `build` folder to
# configure the project before building
option(USE_OPENMP "Enable OpenMP" ON)

# Find and enable OpenMP (not required, if requiring OpenMP,
# use `find_package(OpenMP REQUIRED)` instead)

if (USE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# If build type is "Profile", adding gprof flag (Mac OS does not support gprof)
# The build type can be set be using `cmake -DCMAKE_BUILD_TYPE=Profile ..`
# command in build directory
if (CMAKE_BUILD_TYPE STREQUAL "Profile")
    message(STATUS "BUILD type is Profile. Adding gprof flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

# Add include directory to head file searching path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Gather all source files in src (including subfolders)
file(GLOB SRC_FILES
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*/*.cpp)

# Define executable
add_executable(simulation ${SRC_FILES})
# Link OpenMP if it is found
if (USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(simulation PUBLIC OpenMP::OpenMP_CXX)
endif()

# generate start_simulation.sh file to build directory
configure_file(${PROJECT_SOURCE_DIR}/cmake/start_simulation.sh.in
                ${PROJECT_BINARY_DIR}/start_simulation.sh @ONLY)

# Set the generated start_simulation as executable
file(CHMOD ${PROJECT_BINARY_DIR}/start_simulation.sh
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Define a list of additional clean files
set(ADDITIONAL_CLEAN_FILES
    "${PROJECT_BINARY_DIR}/start_simulation.sh"
    "${PROJECT_BINARY_DIR}/cmake_install.cmake"
    "${PROJECT_BINARY_DIR}/CMakeCache.txt"
    "${PROJECT_BINARY_DIR}/Makefile"
    "${PROJECT_BINARY_DIR}/CMakeFiles"
    "${PROJECT_BINARY_DIR}/compile_commands.json"
)

# Set direcroty property, delete these additional files
# when executing "make clean"
set_directory_properties(PROPERTIES
    ADDITIONAL_CLEAN_FILES "${ADDITIONAL_CLEAN_FILES}"
)


