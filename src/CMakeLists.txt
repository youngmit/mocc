include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mocc-core)

add_subdirectory(mocc-core)

add_executable(mocc mocc.cpp)
target_link_libraries(mocc mocc-core)
