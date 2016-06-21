include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB mc_src "*.cpp")

add_library(mc ${mc_src})
target_link_libraries(mc core pugixml ${HDF5_LIBRARIES} )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
