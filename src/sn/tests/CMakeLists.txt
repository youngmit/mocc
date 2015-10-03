add_executable(testSnBoundary test_sn_boundary.cpp)
target_link_libraries(testSnBoundary mocc-core sn)
add_test(testSnBoundary testSnBoundary)
