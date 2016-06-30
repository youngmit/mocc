if(${BUILD_TESTS})
    add_unit_test(test_teestream)
    add_unit_test(test_fp_utils)
    add_unit_test(test_StringUtils util)
    add_unit_test(test_RNG_LCG)

endif()
