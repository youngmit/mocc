if(${BUILD_TESTS})
add_executable(testAngularQuadrature test_angular_quadrature.cpp)
target_link_libraries(testAngularQuadrature core pugixml)
add_test(testAngularQuadratrue testAngularQuadrature)

add_executable(testMaterial test_material.cpp)
target_link_libraries(testMaterial core pugixml)
add_test(testMaterial testMaterial)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testMaterial)

add_executable(testMesh test_mesh.cpp)
target_link_libraries(testMesh core pugixml)
add_test(testMesh testMesh)

add_executable(testExp test_exponential.cpp)
target_link_libraries(testExp core pugixml)
add_test(testExp testExp)

endif()