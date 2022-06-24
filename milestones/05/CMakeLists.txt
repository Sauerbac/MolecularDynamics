set(milestone05_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/berendsen.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LatticeGen.h)

set(milestone05_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/berendsen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LatticeGen.cpp)

include_directories(../03 ../04)

add_library(milestone05 ${milestone05_HDRS} ${milestone05_SRCS})
target_link_libraries(milestone05 PUBLIC Eigen3::Eigen milestone03 milestone04
                                         ${MPI_CXX_LIBRARIES})

add_executable(milestone05_exe ${milestone05_SRCS} ${milestone05_HDRS})
target_link_libraries(milestone05_exe PUBLIC Eigen3::Eigen milestone03
                                             milestone04 ${MPI_CXX_LIBRARIES})
