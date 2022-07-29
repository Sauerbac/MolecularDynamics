set(milestone07_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/gupta.h ${CMAKE_CURRENT_SOURCE_DIR}/metrics.h)

set(milestone07_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/gupta.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metrics.cpp)

include_directories(../03 ../04 ../05 ../06)

add_library(milestone07 ${milestone07_HDRS} ${milestone07_SRCS})
target_link_libraries(
  milestone07 PUBLIC Eigen3::Eigen milestone03 milestone04 milestone05
                     milestone06 ${MPI_CXX_LIBRARIES})

add_executable(milestone07_exe ${milestone07_HDRS} ${milestone07_SRCS})
target_link_libraries(
  milestone07_exe PUBLIC Eigen3::Eigen milestone03 milestone04 milestone05
                         milestone06 ${MPI_CXX_LIBRARIES})

add_executable(lattice_gen ${CMAKE_CURRENT_SOURCE_DIR}/ih.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/vector.h)
