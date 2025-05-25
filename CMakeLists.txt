
cmake_minimum_required(VERSION 3.10)

# Define project name
project(VorticitySimulation LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directory to head file searching path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Gather all source files in src (including subfolders)
file(GLOB SRC_FILES
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*/*.cpp)

# Define executable
add_executable(simulation ${SRC_FILES})

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
)

# Set direcroty property, delete these additional files when executing "make clean"
set_directory_properties(PROPERTIES
    ADDITIONAL_CLEAN_FILES "${ADDITIONAL_CLEAN_FILES}"
)


