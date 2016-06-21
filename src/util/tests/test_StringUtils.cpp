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

#include <iostream>
#include "util/error.hpp"
#include "util/string_utils.hpp"

TEST(explode_string)
{

    std::string test_string   = " 1 2   3 45 4   \n2 \r7\t 49 ";
    std::vector<int> ref_vec  = {1, 2, 3, 45, 4, 2, 7, 49};
    std::vector<int> test_vec = explode_string<int>(test_string);

    CHECK_EQUAL(ref_vec.size(), test_vec.size());

    for (unsigned i = 0; i < ref_vec.size(); ++i) {
        CHECK_EQUAL(ref_vec[i], test_vec[i]);
    }

    // Test error checking
    test_string = " 1 2 3 4 5 7 3. 1 ";
    CHECK_THROW(explode_string<int>(test_string), mocc::Exception);

    return;
}

TEST(explode_braces)
{
    std::string test_string = "{1 1 1 1 1}{2 2 2  } { 3 3 3 } ";
    auto test_vec = explode_braces(test_string);

    CHECK_EQUAL(3, test_vec.size());

    std::vector<int> tv = {1, 1, 1, 1, 1};
    CHECK_ARRAY_EQUAL(tv, test_vec[0], 5);
    tv = {2, 2, 2};
    CHECK_ARRAY_EQUAL(tv, test_vec[1], 3);
    tv = {3, 3, 3};
    CHECK_ARRAY_EQUAL(tv, test_vec[2], 3);

    test_string = "{ 34 5 2";
    CHECK_THROW(explode_braces(test_string), mocc::Exception);

    test_string = "1 {3 5 1}1 2 3 { 34 5 2} 3";
    test_vec = explode_braces(test_string);
    for(const auto &l: test_vec) {
        for(const auto &v: l) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
    CHECK_EQUAL(1, test_vec[0].size());
    CHECK_EQUAL(1, test_vec[0][0]);
    CHECK_EQUAL(3, test_vec[1].size());
    CHECK_EQUAL(5, test_vec[1][1]);
    CHECK_EQUAL(1, test_vec[2].size());
    CHECK_EQUAL(1, test_vec[3].size());


    test_string   = " 1 2   3 45 4   \n2 \r7\t 49 ";
    test_vec = explode_braces(test_string);
    std::vector<int> flat_vec;
    for(const auto &block: test_vec) {
        for(const auto &v: block){
            flat_vec.push_back(v);
        }
    }
    std::vector<int> ref_vec  = {1, 2, 3, 45, 4, 2, 7, 49};
    CHECK_ARRAY_EQUAL(ref_vec, flat_vec, 8);

    return;
}

int main()
{
    return UnitTest::RunAllTests();
}
