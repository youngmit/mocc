include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/sweepers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/solvers)

add_subdirectory(core)
add_subdirectory(sweepers)
add_subdirectory(solvers)
add_subdirectory(auxiliary)
#add_subdirectory(tests)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/git_SHA1.cpp.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/git_SHA1.cpp" @ONLY)

file(GLOB mocc_src "*.cpp")

add_library(driver "driver.cpp" "git_SHA1.cpp" "input_proc.cpp")

add_executable(mocc "mocc.cpp")
target_link_libraries(mocc driver auxiliary core sweepers solvers ${HDF5_LIBRARIES} ${Blitz_LIBRARY})
