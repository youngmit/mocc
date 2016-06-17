file(GLOB cmdo_src "*.cpp")

add_library(cmdo ${cmdo_src})
target_link_libraries(cmdo core pugixml ${HDF5_LIBRARIES} ${Blitz_LIBRARY})
