if(${BUILD_TESTS})
    add_unit_test(test_Line ${HDF5_LIBRARIES} core)
    add_unit_test(test_Circle ${HDF5_LIBRARIES} core)
    add_unit_test(test_Angle ${HDF5_LIBRARIES} core)

endif()
