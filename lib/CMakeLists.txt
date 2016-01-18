add_library(pugixml pugixml/src/pugixml.cpp)

include_directories("${CMAKE_CURRENT_SRC_DIRECTORY}/pugixml/include"
    "${CMAKE_CURRENT_SRC_DIRECTORY}/")

add_subdirectory(unittest-cpp)
