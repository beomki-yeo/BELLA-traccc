# Set up the project.
cmake_minimum_required( VERSION 3.22 )
project( momentum_resolution LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 20 CACHE STRING "The (host) C++ standard to use" )

find_package( Boost 1.86.0 REQUIRED COMPONENTS program_options filesystem log)

# Include covfie and traccc
add_subdirectory(extern/covfie)
add_subdirectory(extern/traccc)

# Write B field
add_executable( write_bfield src/write_bfield.cpp )
target_link_libraries( write_bfield )

# Build simulation
add_executable( do_telescope_simulation src/telescope_simulation.cpp )
target_link_libraries( do_telescope_simulation PRIVATE 
    vecmem::core detray::io detray::detectors 
    traccc::core traccc::io traccc::options traccc::simulation 
    Boost::filesystem )

# Build fitting example
add_executable( do_truth_fitting_momentum_residual src/truth_fitting_momentum_residual.cpp )
target_link_libraries( do_truth_fitting_momentum_residual PRIVATE 
    vecmem::core detray::io detray::detectors 
    traccc::core traccc::io traccc::options traccc::performance )