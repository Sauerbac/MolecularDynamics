set(milestone09_HDRS # Add header files here
    )

set(milestone09_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

include_directories(../03 ../04 ../05 ../06 ../07 ../08)

add_library(milestone09 ${milestone09_HDRS} ${milestone09_SRCS})
target_link_libraries(
  milestone09
  PUBLIC Eigen3::Eigen
         milestone03
         milestone04
         milestone05
         milestone06
         milestone07
         milestone08
         ${MPI_CXX_LIBRARIES})

add_executable(milestone09_exe ${milestone09_HDRS} ${milestone09_SRCS})
target_link_libraries(
  milestone09_exe
  PUBLIC Eigen3::Eigen
         milestone03
         milestone04
         milestone05
         milestone06
         milestone07
         milestone08
         ${MPI_CXX_LIBRARIES})
