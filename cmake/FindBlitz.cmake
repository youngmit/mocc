#
# Locate Blitz++++ include paths and libraries
# Blitz++++ can be found at http://www.oonumerics.org/Blitz++/
# Written by Michael Hammer, michael.hammer_at_tugraz.at

# This module defines
# Blitz++_INCLUDE_DIR, where to find ptlib.h, etc.
# Blitz++_LIBRARIES, the libraries to link against to use pwlib.
# Blitz++_FOUND, If false, don't try to use pwlib.

include(LibFindMacros)

FIND_PATH(Blitz_INCLUDE_DIR blitz/blitz.h
    "$ENV{Blitz}/include"
    /usr/local/include
    /usr/include
)

FIND_LIBRARY(Blitz_LIBRARY libblitz.a
  PATHS
    "$ENV{Blitz}/lib"
    /usr/local/lib
    /usr/lib
)

MESSAGE(STATUS "Blitz env: $ENV{Blitz}")

include(FindPackageHandleStandardArgs)
MESSAGE(STATUS "Looking for Blitz++")
find_package_handle_standard_args(Blitz DEFAULT_MSG
        Blitz_LIBRARY Blitz_INCLUDE_DIR)
MESSAGE(STATUS "Found Blitz++")
MESSAGE(STATUS ${Blitz_INCLUDE_DIR})
MESSAGE(STATUS ${Blitz_LIBRARY})

MARK_AS_ADVANCED(
  Blitz_INCLUDE_DIR
  Blitz_LIBRARIES
) 
