include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB mocc_src "*.cpp")

add_executable(mocc ${mocc_src} ${pugixml})
target_link_libraries(mocc pugixml)
target_link_libraries(mocc ${Boost_LIBRARIES})
