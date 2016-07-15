include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(sn)
add_subdirectory(moc)
add_subdirectory(cmdo)
add_subdirectory(mc)

file(GLOB sweepers_src "*.cpp")

add_library(sweepers ${sweepers_src})
target_link_libraries(sweepers core sn moc cmdo mc pugixml )
