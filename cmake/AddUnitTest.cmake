#
# add_unit_test( <target> <link targets>... )
#
# Add a unit test, which should use the UnitTest++ harness. The first argument
# should be the name of the unit test, which must be the same as the source
# file, without its '.cpp' extension. Following arguments are the libraries
# against which the test should be linked.
#
function(add_unit_test target)
        add_executable(${target} ${target}.cpp )
    target_link_libraries(${target} ${ARGN} UnitTest++)

    add_test(${target} ${target})
endfunction()
