if(${BUILD_TESTS})
    add_unit_test(test_h5file ${HDF5_LIBRARIES} core)

    add_unit_test(test_AngularQuadrature ${HDF5_LIBRARIES} core)

    add_unit_test(test_Assembly ${HDF5_LIBRARIES} core)

    add_unit_test(test_PinMesh_Cyl ${HDF5_LIBRARIES} core)

    add_unit_test(test_PinMesh_Rect ${HDF5_LIBRARIES} core)

    add_unit_test(test_ScatteringMatrix ${HDF5_LIBRARIES} core)

    add_unit_test(test_Material core pugixml)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl
        ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_Material)

    add_unit_test(test_Mesh core pugixml)
    add_unit_test(test_CoreMesh core pugixml)

    add_unit_test(test_Exponential core pugixml)

    add_unit_test(test_BoundaryCondition core pugixml ${HDF5_LIBRARIES})

    add_unit_test(test_CMFD core pugixml)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/3x5.xml
        ${CMAKE_CURRENT_BINARY_DIR}/3x5.xml test_CMFD)

    add_unit_test(test_XSMeshHomogenized core pugixml ${HDF5_LIBRARIES})
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl
            ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_XSMeshHomogenized)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/2x3_stack.xml
            ${CMAKE_CURRENT_BINARY_DIR}/2x3_stack.xml test_XSMeshHomogenized)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/2x3_1.xml
            ${CMAKE_CURRENT_BINARY_DIR}/2x3_1.xml test_XSMeshHomogenized)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/2x3_2.xml
            ${CMAKE_CURRENT_BINARY_DIR}/2x3_2.xml test_XSMeshHomogenized)

    add_unit_test(test_XSMesh core pugixml ${HDF5_LIBRARIES})

endif()
