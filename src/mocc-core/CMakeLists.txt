include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

file(GLOB core_src "*.cpp")

add_library(mocc-core ${core_src})
target_link_libraries(mocc-core pugixml hdf5_cpp)
target_link_libraries(mocc-core ${Boost_LIBRARIES})
