include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE)

set(MolecularDynamics_tests_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/test_hello_world.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_neighbors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_lj_direct_summation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_verlet.cpp)

# Tests
include_directories(
  ${CMAKE_SOURCE_DIR}/milestones/03 ${CMAKE_SOURCE_DIR}/milestones/04
  ${CMAKE_SOURCE_DIR}/milestones/05 ${CMAKE_SOURCE_DIR}/milestones/06)
add_executable(MolecularDynamics_tests ${MolecularDynamics_tests_SRCS})

target_link_libraries(
  MolecularDynamics_tests
  PUBLIC Eigen3::Eigen
         gtest
         gtest_main
         milestone03
         milestone04
         milestone05
         milestone06)

include(GoogleTest)
gtest_discover_tests(MolecularDynamics_tests)

add_test(NAME MolecularDynamics_tests COMMAND MolecularDynamics_tests)
