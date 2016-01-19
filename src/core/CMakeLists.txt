include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

file(GLOB core_src "*.cpp")

add_library(core ${core_src})
target_link_libraries(core pugixml hdf5 )
target_link_libraries(core)
