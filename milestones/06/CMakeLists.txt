set(milestone06_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/neighbors.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ForcesEnergiesCutoff.h)

set(milestone06_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/ForcesEnergiesCutoff.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/neighbors.cpp)

include_directories(../03 ../04 ../05)

add_library(milestone06 ${milestone06_HDRS} ${milestone06_SRCS})
target_link_libraries(milestone06 PUBLIC Eigen3::Eigen milestone03 milestone04
                                         milestone05 ${MPI_CXX_LIBRARIES})

add_executable(milestone06_exe ${milestone06_HDRS} ${milestone06_SRCS})
target_link_libraries(
  milestone06_exe PUBLIC Eigen3::Eigen milestone03 milestone04 milestone05
                         ${MPI_CXX_LIBRARIES})
