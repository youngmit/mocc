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

#pragma once

#include <memory>

#include "core/output_interface.hpp"
#include "core/transport_sweeper.hpp"

namespace mocc {
/**
* This provides a virtual base type, which shall provide a solve() and
* step() method. At the highest level of the heirarchy, the driver calls
* solve() and that should invoke everything that is necessary.
*/
class Solver : public HasOutput {
public:
    virtual ~Solver(){};

    /**
    * Perform a full solution to the class of problem that the most-derived
    * Solver type is designed to solve. This is usually called upon the top-
    * level Solver by the driver.
    */
    virtual void solve() = 0;

    /**
    * Perfom some sort of intermediate step in solving the problem of
    * interest, typically as part of another solver. What specifically is
    * done is quite solver specific, so check the derived class to see what
    * it does for a specific case.
    */
    virtual void step() = 0;

    /**
    * Return a pointer to a transport sweeper object. If the Solver does
    * not actually have a sweeper, return nullptr.
    */
    virtual const TransportSweeper *sweeper() const
    {
        return nullptr;
    }

private:
};

typedef std::shared_ptr<Solver> SP_Solver_t;
}
