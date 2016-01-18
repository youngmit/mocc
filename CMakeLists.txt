project(MOCC)
cmake_minimum_required(VERSION 3.2)

SET(BUILD_TESTS true CACHE BOOL "Enable compilation of tests")
SET(PROFILE false CACHE BOOL "Enable profiling")
MESSAGE(STATUS "Tests: ${BUILD_TESTS}")
MESSAGE(STATUS "Profiling: ${PROFILE}")
enable_testing()

MACRO(COPY_FILE_IF_CHANGED in_file out_file target)
    message("Copying file: ${in_file} to: ${out_file}")
    IF(${in_file} IS_NEWER_THAN ${out_file})    
        ADD_CUSTOM_COMMAND (
            TARGET     ${target}
            POST_BUILD
            COMMAND    ${CMAKE_COMMAND}
            ARGS       -E copy ${in_file} ${out_file}
        )
    ENDIF(${in_file} IS_NEWER_THAN ${out_file})
ENDMACRO(COPY_FILE_IF_CHANGED)


message("Compiler ID: '${CMAKE_CXX_COMPILER_ID}'")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Using GNU")
    set(warnings "-Wall")
    set(options "-std=c++14 -fopenmp")
    if (${PROFILE})
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
        SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(warnings "/W4 /XW /EHsc")
else()
    set(warnings "-Wall")
    set(options "-std=c++14")
endif()

message("Flags: ${warnings} ${options}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings} ${options}"
    CACHE STRING "Flags used by the compiler during all build types." FORCE)
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings} ${options}"
    CACHE STRING "Flags used by the compiler during all build types." FORCE)
#set(CMAKE_CXX_LINK_FLAGS ${link})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(HDF5)

find_package(Blitz REQUIRED)

include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/pugixml/src")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/unittest-cpp")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")
include_directories(SYSTEM ${Blitz_INCLUDE_DIR})


include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
add_subdirectory(lib)

add_subdirectory(src)
