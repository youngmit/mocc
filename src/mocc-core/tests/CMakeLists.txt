if(${BUILD_TESTS})
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

endif()
