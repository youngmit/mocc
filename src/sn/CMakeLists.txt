file(GLOB sn_src "*.cpp")

add_library(sn ${sn_src})
target_link_libraries(sn mocc-core pugixml hdf5_cpp)
