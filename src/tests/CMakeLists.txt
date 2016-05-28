if(${BUILD_TESTS})
        set(link_tests "driver;auxiliary;core;sweepers;solvers;${HDF5_LIBRARIES};pugixml")

        add_unit_test(test_3x3 ${link_tests})
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/3x3.xml
            ${CMAKE_CURRENT_BINARY_DIR}/3x3.xml test_3x3)
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl
            ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_3x3)




endif()
