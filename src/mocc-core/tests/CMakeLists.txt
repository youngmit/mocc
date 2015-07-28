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
