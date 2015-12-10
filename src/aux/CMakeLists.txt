file(GLOB aux_src "*.cpp")

add_library(aux ${aux_src})
target_link_libraries(aux pugixml ${Blitz_LIBRARY})
