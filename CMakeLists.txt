project(MOCC)
cmake_minimum_required(VERSION 3.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

SET(BUILD_TESTS true CACHE BOOL "Enable compilation of tests")
MESSAGE(STATUS "Tests: ${BUILD_TESTS}")
SET(PROFILE false CACHE BOOL "Enable profiling")
MESSAGE(STATUS "Profiling: ${PROFILE}")
SET(COVERAGE false CACHE BOOL "Enable code coverage instrumentation")
MESSAGE(STATUS "Coverage: ${COVERAGE}")

enable_testing()

include(AddUnitTest)

CONFIGURE_FILE("${CMAKE_MODULE_PATH}/CTestCustom.cmake" "${CMAKE_BINARY_DIR}/CTestCustom.cmake")

include(CTest)

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
        MESSAGE(STATUS "Release flags were: ${CMAKE_CXX_FLAGS_RELEASE}")
        SET(CMAKE_CXX_FLAGS_RELEASE " -O1 -DNDEBUG")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
        SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
    endif()
    if (${COVERAGE})
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(warnings "/W4 /XW /EHsc")
else()
    set(warnings "-Wall")
    set(options "-std=c++14")
endif()

message("Flags: ${warnings} ${options}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings} ${options}" )
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings} ${options}")

# What a shame...
if(WIN32)
    # windows cmake internal lookups are broken for now
    # will lookup for headers and shared libs given HDF_DIR env
    find_path(HDF5_INCLUDE_DIRS hdf5.h HINTS "$ENV{HDF5_DIR}\\..\\include")
    find_library(HDF5_ZLIB_LIBRARY NAMES libzlib HINTS "$ENV{HDF5_DIR}\\..\\lib")
    find_library(HDF5_SZIP_LIBRARY NAMES libszip HINTS "$ENV{HDF5_DIR}\\..\\lib")
    find_library(HDF5_C_LIBRARY NAMES libhdf5 HINTS "$ENV{HDF5_DIR}\\..\\lib")
    find_library(HDF5_CXX_LIBRARY NAMES libhdf5_cpp HINTS "$ENV{HDF5_DIR}\\..\\lib")
    if(HDF5_INCLUDE_DIRS AND HDF5_C_LIBRARY)
        MESSAGE(STATUS "Found HDF5!")
        set(HDF5_FOUND "YES")
        set(HDF5_LIBRARIES ${HDF5_SZIP_LIBRARY}
                           ${HDF5_ZLIB_LIBRARY}
                           ${HDF5_C_LIBRARY}
                           ${HDF5_CXX_LIBRARY})
        mark_as_advanced(HDF5_LIBRARIES)
        include_directories(${HDF5_INCLUDE_DIRS})
    else()
        message(STATUS "Could not find HDF5")
        set(HDF5_FOUND "NO")
    endif()
    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_C_LIBRARY: ${HDF5_C_LIBRARY}")
    message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
else()
    find_package(HDF5)
endif()

find_package(Blitz REQUIRED)

include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/pugixml/src")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/unittest-cpp")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")
include_directories(SYSTEM ${Blitz_INCLUDE_DIR})

add_subdirectory(lib)

add_subdirectory(src)
