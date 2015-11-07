#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <iostream>

#include "mocc-core/exponential.hpp"

using namespace mocc;
using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE( testall )
{

    Exponential_Linear exp;
    BOOST_CHECK_EQUAL( exp.exp(0.0), 1.0 );
    BOOST_CHECK_EQUAL( exp.exp(-1.0), 0.36787944117144233 );
    BOOST_CHECK_EQUAL( exp.exp(-0.314), 0.7305190281594249 );
    cout << "Max error: " << exp.max_error() << endl;
}
