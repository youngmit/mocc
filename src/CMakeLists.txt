include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/sweepers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/solvers)

add_subdirectory(core)
add_subdirectory(sweepers)
add_subdirectory(solvers)
add_subdirectory(auxiliary)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/git_SHA1.cpp.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/git_SHA1.cpp" @ONLY)

file(GLOB mocc_src "*.cpp")

add_executable(mocc ${mocc_src})
target_link_libraries(mocc auxiliary core sweepers solvers hdf5_cpp ${Blitz_LIBRARY})
