#undef __STRICT_ANSI__
#undef _REENT_ONLY
#define BOOST_TEST_MAIN
#include <stdlib.h>
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <iostream>

#include "core/exponential.hpp"

using namespace mocc;
using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE( testall )
{

    Exponential_Linear exp;
    cout << "Max error: " << exp.max_error() << endl;

    real_t max_err = 0.0;
    for( real_t x=-10.0; x<0.0; x+=0.1 ) {
        real_t exp_t = exp.exp(x);
        real_t exp_r = std::exp(x);
        max_err = std::max(max_err, std::abs(exp_r - exp_t));
        cout << x << " " << exp_r << " " << exp_t << endl;
        BOOST_CHECK( std::abs(exp_r - exp_t) < 1e-9 );
    }

    cout << "max_err: " << max_err << endl;
}
