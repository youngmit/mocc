if(${BUILD_TESTS})
    add_executable(testH5File test_h5file.cpp)
    target_link_libraries(testH5File core hdf5_cpp UnitTest++)
    add_test(testH5File testH5File)

    add_executable(testAngularQuadrature test_angular_quadrature.cpp)
    target_link_libraries(testAngularQuadrature core pugixml UnitTest++)
    add_test(testAngularQuadratrue testAngularQuadrature)
    
    add_executable(testMaterial test_material.cpp)
    target_link_libraries(testMaterial core pugixml UnitTest++)
    add_test(testMaterial testMaterial)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
        ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testMaterial)
    
    add_executable(testMesh test_mesh.cpp)
    target_link_libraries(testMesh core pugixml UnitTest++)
    add_test(testMesh testMesh)
    
    add_executable(testExp test_exponential.cpp)
    target_link_libraries(testExp core pugixml UnitTest++)
    add_test(testExp testExp)
    
    add_executable(testXSMeshHom test_xsmeshhomogenized.cpp)
    target_link_libraries(testXSMeshHom core pugixml hdf5_cpp UnitTest++)
    add_test(testXSMeshHom testXSMeshHom)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
            ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testXSMeshHom)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
            ${CMAKE_CURRENT_BINARY_DIR}/square.xml testXSMeshHom)
    copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square_2.xml 
            ${CMAKE_CURRENT_BINARY_DIR}/square_2.xml testXSMeshHom)
endif()
