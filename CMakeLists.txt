cmake_minimum_required(VERSION 3.10)

# Project name
project(MolFFSim LANGUAGES CXX)

# Set C++ compiler (replace "g++" with your desired compiler)
set(CMAKE_CXX_COMPILER g++)
    
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Xpreprocessor -fopenmp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_EXE_LINKER_FLAGS "-lomp")
    set(OpenMP_INCLUDE_DIRS "/opt/homebrew/opt/libomp/include")
    set(OpenMP_LIBRARIES "/opt/homebrew/opt/libomp/lib/libomp.dylib")
else()
    # Check if OpenMP is available
    find_package(OpenMP REQUIRED)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif() 

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen library
find_package (Eigen3 3.3 NO_MODULE)

# Find autodiff library
find_package(autodiff REQUIRED)

# Find liblbfgs library
find_package(LBFGS REQUIRED)

# Add /usr/local/include to the list of include directories
include_directories(/usr/local/include)

# Include eigen3 headers
include_directories(${EIGEN3_INCLUDE_DIRS})

# Include autodiff headers
include_directories(${AUTODIFF_INCLUDE_DIRS})

# Include lbfgs headers
include_directories(${LBFGS_INCLUDE_DIRS})

# Add an executable target
file(GLOB SOURCES "src/*.cpp" "src/*.hpp")
add_executable(MolFFSim ${SOURCES})

# Add link with OpenMP
if(OpenMP_FOUND)
    target_link_libraries(MolFFSim PRIVATE OpenMP::OpenMP_CXX)
else()
    target_include_directories(MolFFSim PRIVATE ${OpenMP_INCLUDE_DIRS})
    target_link_libraries(MolFFSim PRIVATE ${OpenMP_LIBRARIES})
endif()

# Link the eigen library
target_link_libraries(MolFFSim ${EIGEN3_LIBRARIES})

# Link the autodiff library
target_link_libraries(MolFFSim PRIVATE ${AUTODIFF_LIBRARIES})

# Link the liblbfgs library
target_link_libraries(MolFFSim PRIVATE ${LBFGS_LIBRARIES})

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