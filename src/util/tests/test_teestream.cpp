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

#include "util/tee_stream.hpp"

#include <iostream>
#include <sstream>

TEST(test_TeeStream)
{

    std::stringstream truth;

    std::stringstream s1;
    std::stringstream s2;

    // Make sure the output is duplicated across streams
    TeeStream ts(s1, s2);

    ts << "something about foxes" << std::endl;
    truth << "something about foxes" << std::endl;

    CHECK_EQUAL(truth.str(), s1.str());
    CHECK_EQUAL(truth.str(), s2.str());

    std::cout << s1.str();
    std::cout << s2.str();

    std::stringstream s3;
    std::stringstream s4;

    // Make sure the new streams are used
    ts.reset(s3, s4);

    ts << "something else about a lazy dog" << std::endl;
    truth.str("");
    truth << "something else about a lazy dog" << std::endl;

    std::cout << s3.str();
    std::cout << s4.str();

    CHECK_EQUAL(truth.str(), s3.str());
    CHECK_EQUAL(truth.str(), s4.str());

    onullstream null_stream;

    s1.str("");

    // Make sure the onullstream works right
    ts.reset(s1, null_stream);

    ts << "such jump. very wow." << std::endl;
    truth.str("");
    truth << "such jump. very wow." << std::endl;
    CHECK_EQUAL(s1.str(), truth.str());

    // Make sure that the previously-associated stream (s3) is unaltered.
    CHECK_EQUAL(s3.str(), s4.str());

    std::cout << s1.str();
}

int main()
{
    return UnitTest::RunAllTests();
}
