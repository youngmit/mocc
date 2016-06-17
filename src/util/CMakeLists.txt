include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Commit the force inline stuff to the force_inline.hpp header so that it can be
# included where needed and used to do a platform-independent forced inline
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/force_inline.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/force_inline.hpp")



add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
