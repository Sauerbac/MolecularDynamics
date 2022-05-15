cmake_minimum_required(VERSION 3.11)

project(MolecularDynamics)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(MolecularDynamics_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/Atoms.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Types.h
        )

set(MolecularDynamics_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/Atoms.cpp
        )

set(MolecularDynamics_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )

add_subdirectory(tests)

# Simulation code
add_executable(MolecularDynamics ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${MolecularDynamics_executable_SRCS})
target_link_libraries(MolecularDynamics PUBLIC Eigen3::Eigen)
