project(MOCC)
cmake_minimum_required(VERSION 3.2)

message("Compiler ID: '${CMAKE_CXX_COMPILER_ID}'")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Using GNU")
    set(warnings "-Wall")
    set(options "-std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(warnings "/W4 /XW /EHsc")
else ()
    message("wtf")
endif()
    set(warnings "-Wall")
    set(options "-std=c++11")
message("Compiler ID: '${CMAKE_CXX_COMPILER_ID}'")

message("Flags: ${warnings} ${options}")

#if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${warnings} ${options}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${warnings} ${options}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
#endif()

find_package(Boost REQUIRED COMPONENTS regex)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib/pugixml/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
add_subdirectory(lib)

add_subdirectory(src)
