if(${BUILD_TESTS})
add_executable(testRayData test_ray_data.cpp)
target_link_libraries(testRayData mocc-core pugixml)
add_test(testRayData testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/square.xml testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testRayData)

add_executable(testAngularQuadrature test_angular_quadrature.cpp)
target_link_libraries(testAngularQuadrature mocc-core pugixml)
add_test(testAngularQuadratrue testAngularQuadrature)

add_executable(testMaterial test_material.cpp)
target_link_libraries(testMaterial mocc-core pugixml)
add_test(testMaterial testMaterial)

add_executable(testMesh test_mesh.cpp)
target_link_libraries(testMesh mocc-core pugixml)
add_test(testMesh testMesh)

add_executable(testExp test_exponential.cpp)
target_link_libraries(testExp mocc-core pugixml)
add_test(testExp testExp)


add_executable(testRay test_ray.cpp)
target_link_libraries(testRay mocc-core pugixml)
add_test(testRay testRay)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/6x5.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/6x5.xml testRay)
endif()