include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

file(GLOB core_geom_src "*.cpp")

message(${core_geom_src})

set( core_src ${core_src} ${core_geom_src} PARENT_SCOPE)

