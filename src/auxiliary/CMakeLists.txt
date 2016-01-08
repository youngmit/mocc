file(GLOB aux_src "*.cpp")

add_library(auxiliary ${aux_src})
target_link_libraries(auxiliary pugixml ${Blitz_LIBRARY})
