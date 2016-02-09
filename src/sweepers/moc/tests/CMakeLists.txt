if(${BUILD_TESTS})

        
add_executable(testRayData test_ray_data.cpp)
target_link_libraries(testRayData core moc pugixml UnitTest++)
add_test(testRayData testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/square.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/square.xml testRayData)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl 
    ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl testRayData)

add_executable(testRay test_ray.cpp)
target_link_libraries(testRay core moc pugixml UnitTest++)
add_test(testRay testRay)
copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/6x5.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/6x5.xml testRay)


add_executable(testMoC_IHM test_moc_ihm.cpp)
target_link_libraries(testMoC_IHM core moc pugixml UnitTest++)
add_test(testMoC_IHM testMoC_IHM)
endif()
