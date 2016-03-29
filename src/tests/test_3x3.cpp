#include "UnitTest++/UnitTest++.h"

#include "driver.hpp"

TEST(test_3x3) {
    run("3x3.xml");
}

int main() {
    return UnitTest::RunAllTests();
}
