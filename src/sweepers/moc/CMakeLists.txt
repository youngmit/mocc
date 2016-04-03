include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

file(GLOB moc_src "*.cpp")

add_library(moc ${moc_src})
target_link_libraries(moc core pugixml ${HDF5_LIBRARIES} )
target_link_libraries(moc)
