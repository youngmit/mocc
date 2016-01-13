include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB solvers_src "*.cpp")

add_library(solvers ${solvers_src})
target_link_libraries(solvers sweepers pugixml hdf5_cpp )
