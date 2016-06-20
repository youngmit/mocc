if(${BUILD_TESTS})

    add_unit_test(test_ParticlePusher core mc pugixml)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/1g.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/1g.xsl test_ParticlePusher)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/tunnel.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/tunnel.xml test_ParticlePusher)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/square.xml test_ParticlePusher)
endif()
