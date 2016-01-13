include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(sn)
add_subdirectory(moc)
add_subdirectory(cmdo)


file(GLOB sweepers_src "*.cpp")

add_library(sweepers ${sweepers_src})
target_link_libraries(sweepers mocc-core sn moc cmdo pugixml )
