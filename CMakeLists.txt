cmake_minimum_required(VERSION 3.10)

# Project name
project(MolFFSim LANGUAGES CXX)

# Set C++ compiler (replace "g++" with your desired compiler)
set(CMAKE_CXX_COMPILER g++)

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




