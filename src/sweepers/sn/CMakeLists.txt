file(GLOB sn_src "*.cpp")

add_library(sn ${sn_src})
target_link_libraries(sn core pugixml ${HDF5_LIBRARIES})
