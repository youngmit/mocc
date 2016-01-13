if(${BUILD_TESTS})

        
add_executable(testRayData test_ray_data.cpp)
target_link_libraries(testRayData mocc-core moc pugixml)
add_test(testRayData testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/square.xml testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testRayData)

add_executable(testRay test_ray.cpp)
target_link_libraries(testRay mocc-core moc pugixml)
add_test(testRay testRay)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/6x5.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/6x5.xml testRay)
endif()
