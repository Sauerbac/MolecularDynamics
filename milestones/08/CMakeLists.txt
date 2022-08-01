set(milestone08_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi_support.h
    ${CMAKE_CURRENT_SOURCE_DIR}/domain.h)

set(milestone08_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/domain.cpp)

include_directories(../03 ../04 ../05 ../06 ../07)

add_library(milestone08 ${milestone08_HDRS} ${milestone08_SRCS})
target_link_libraries(
  milestone08
  PUBLIC Eigen3::Eigen
         milestone03
         milestone04
         milestone05
         milestone06
         milestone07
         ${MPI_CXX_LIBRARIES})

add_executable(milestone08_exe ${milestone08_HDRS} ${milestone08_SRCS})
target_link_libraries(
  milestone08_exe
  PUBLIC Eigen3::Eigen
         milestone03
         milestone04
         milestone05
         milestone06
         milestone07
         ${MPI_CXX_LIBRARIES})
