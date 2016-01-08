include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mocc-core)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/moc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/sn)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cmdo)

add_subdirectory(mocc-core)
add_subdirectory(moc)
add_subdirectory(sn)
add_subdirectory(cmdo)
add_subdirectory(auxiliary)

file(GLOB mocc_src "*.cpp")

add_executable(mocc ${mocc_src})
target_link_libraries(mocc auxiliary mocc-core moc sn cmdo hdf5 ${Blitz_LIBRARY})
