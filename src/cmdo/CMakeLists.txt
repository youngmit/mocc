file(GLOB cmdo_src "*.cpp")

add_library(cmdo ${cmdo_src})
target_link_libraries(cmdo pugixml hdf5_cpp ${Blitz_LIBRARY})
