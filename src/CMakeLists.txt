include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mocc-core)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cmdo)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/sn)

add_subdirectory(mocc-core)
add_subdirectory(sn)
add_subdirectory(cmdo)

file(GLOB mocc_src "*.cpp")

add_executable(mocc ${mocc_src})
target_link_libraries(mocc mocc-core sn cmdo hdf5 ${Blitz_LIBRARY})
