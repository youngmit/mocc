/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "UnitTest++/UnitTest++.h"

#include "blitz_typedefs.hpp"
#include "error.hpp"
#include "h5file.hpp"

using std::cout;
using std::endl;

using namespace mocc;

TEST(test_write)
{
    // Make an HDF5 file
    H5Node h5test("test_file.h5", H5Access::WRITE);

    ArrayB3 data(3, 7, 14);
    data = 0.0;
    data(0, 0, 0) = 1.0;
    data(2, 4, 8) = 7.3;

    h5test.write("test_array", data);

    H5Node g = h5test.create_group("group_a");

    data(2, 2, 6) = 11.32;
    g.write("d", data);

    ArrayB1 d1(10);
    d1    = 1.0;
    d1(3) = 5.32;
    g.write("one_d", d1);

    g = h5test.create_group("/group_b");
    g = h5test.create_group("group_a/sub_1");
    g = h5test.create_group("/group_a/sub_2");

    g.create_group("/foo");
}

TEST(test_read)
{
    H5Node h5test("test_file.h5", H5Access::READ);

    ArrayB3 data;
    h5test.read("test_array", data);
    CHECK_EQUAL(data(0, 0, 0), 1.0);
    CHECK_EQUAL(data(2, 4, 8), 7.3);

    ArrayB1 d1(11);
    CHECK_THROW(h5test.read("/group_a/one_d", d1), Exception);
    d1.resize(10);
    h5test.read("/group_a/one_d", d1);
    CHECK_EQUAL(5.32, d1(3));

    // make sure we arent allowing write-like operations
    CHECK_THROW(h5test.create_group("falala"), Exception);
    // CHECK_THROW(h5test.write, Exception);
}

int main(int, const char *[])
{
    return UnitTest::RunAllTests();
}
