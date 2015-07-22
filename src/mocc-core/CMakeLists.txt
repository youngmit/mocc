include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB mocc_src "*.cpp")

add_library(mocc-core ${mocc_src})
target_link_libraries(mocc-core pugixml)
target_link_libraries(mocc-core ${Boost_LIBRARIES})
