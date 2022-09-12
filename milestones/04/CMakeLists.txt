set(milestone04_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/ForcesEnergies.h
    ${CMAKE_CURRENT_SOURCE_DIR}/xyz.h ${CMAKE_CURRENT_SOURCE_DIR}/metrics.h)

set(milestone04_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/xyz.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ForcesEnergies.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metrics.cpp)

include_directories(../03)

add_library(milestone04 ${milestone04_HDRS} ${milestone04_SRCS})
target_link_libraries(milestone04 PUBLIC Eigen3::Eigen milestone03
                                         ${MPI_CXX_LIBRARIES})

add_executable(milestone04_exe ${milestone04_SRCS} ${milestone04_HDRS})
target_link_libraries(milestone04_exe PUBLIC Eigen3::Eigen milestone03
                                             ${MPI_CXX_LIBRARIES})
