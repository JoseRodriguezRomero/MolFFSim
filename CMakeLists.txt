cmake_minimum_required(VERSION 3.10)

# Project name
project(MolFFSim LANGUAGES CXX)

# Set C++ compiler (replace "g++" with your desired compiler)
set(CMAKE_CXX_COMPILER g++-13)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen library
find_package(eigen3 REQUIRED)

# Find autodiff library
find_package(autodiff REQUIRED)

# Add /usr/local/include to the list of include directories
include_directories(/usr/local/include)

# Include eigen3 headers
include_directories(${EIGEN_INCLUDE_DIRS})

# Include directories for autodiff
include_directories(${AUTODIFF_INCLUDE_DIRS})

# Add an executable target
file(GLOB SOURCES "src/*.cpp" "src/*.hpp")
add_executable(MolFFSim ${SOURCES})

# Link the eigen library
target_link_libraries(MolFFSim PRIVATE ${EIGEN_LIBRARIES})

# Link the autodiff library
target_link_libraries(MolFFSim PRIVATE ${AUTODIFF_LIBRARIES})

# Install MolFFSim
install(
    TARGETS MolFFSim
    DESTINATION /usr/local/share/MolFFSim
)

install(CODE "execute_process(COMMAND ln -sf /usr/local/share/MolFFSim/MolFFSim /usr/local/bin/MolFFSim)")

# Install the atomic basis functions data files
file(GLOB atomic_basis "${CMAKE_SOURCE_DIR}/ATOMIC_BASIS/*.txt")
install(
    FILES ${atomic_basis}
    DESTINATION /usr/local/share/MolFFSim/ATOMIC_BASIS/
)

# Install the exchange and correlation coefficient data files
file(GLOB xc_coefficients "${CMAKE_SOURCE_DIR}/XC_COEFFICIENTS/*.txt")
install(
    FILES ${xc_coefficients}
    DESTINATION /usr/local/share/MolFFSim/XC_COEFFICIENTS/
)

