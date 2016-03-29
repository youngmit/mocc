if(${BUILD_TESTS})

add_unit_test(test_RayData core moc pugixml)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/square.xml test_RayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_RayData)

add_unit_test(test_Ray core moc pugixml)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/6x5.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/6x5.xml test_Ray)

add_unit_test(test_MoC_IHM core moc pugixml)
endif()
