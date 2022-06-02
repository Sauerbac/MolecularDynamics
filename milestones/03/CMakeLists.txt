set(milestone03_HDRS # Add header files here
    ${CMAKE_CURRENT_SOURCE_DIR}/Atoms.h ${CMAKE_CURRENT_SOURCE_DIR}/Types.h)

set(milestone03_SRCS # Add source files here
    ${CMAKE_CURRENT_SOURCE_DIR}/Atoms.cpp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(milestone03 ${milestone03_HDRS} ${milestone03_SRCS})
target_link_libraries(milestone03 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})

add_executable(milestone03_exe ${milestone03_SRCS} ${milestone03_HDRS})
target_link_libraries(milestone03_exe PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})
