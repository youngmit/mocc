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

#include <iostream>

#include "pugixml.hpp"

#include "core/eigen_interface.hpp"
#include "core/core_mesh.hpp"
#include "core/h5file.hpp"
#include "core/cmfd.hpp"
#include "core/transport_sweeper.hpp"

#include "fixed_source_solver.hpp"
#include "solver.hpp"

namespace mocc{
    enum class CMFDConvergence {
        FIXED, // Converge the cmfd to a fixed set of convergence criteria
        FLOAT
    };


    struct ConvergenceCriteria {
        ConvergenceCriteria( real_t k, real_t error_k, real_t error_psi ):
            k(k),
            error_k(error_k),
            error_psi(error_psi) { }
        real_t k;
        real_t error_k;
        real_t error_psi;

        // Stream insertion
        friend std::ostream& operator<<( std::ostream& os,
                ConvergenceCriteria conv);
    };


    class EigenSolver: public Solver {
    public:
        EigenSolver( const pugi::xml_node &input, const CoreMesh &mesh );
        void solve();
        void step();

        const TransportSweeper* sweeper() const {
            return fss_.sweeper();
        }

        // Implement the output interface
        void output( H5Node &file ) const;

    private:
        // Data
        FixedSourceSolver fss_;

        // Fission source, and previous iterate
        ArrayB1 fission_source_;
        ArrayB1 fission_source_prev_;

        // Current guess for k
        real_t keff_;

        // Previous guess for k
        real_t keff_prev_;

        // Convergence criterion for the system eigenvalue
        real_t tolerance_k_;

        // Convergence criterion for the fission source distribution (L-2 norm)
        real_t tolerance_psi_;

        real_t error_k_;
        real_t error_psi_;

        // Maximum allowable outer iterations
        unsigned int max_iterations_;
        unsigned int min_iterations_;

        // Vector of the convergence criteria. We will export these to the HDF5
        // file at the end of the run for posterity
        std::vector<ConvergenceCriteria> convergence_;

        // CMFD accelerator
        UP_CMFD_t cmfd_;


        // Methods
        // Print the current state of the eigenvalue solver
        void print( int iter, ConvergenceCriteria conv );

        /**
         * \brief Perform a CMFD accelerator solve
         */
        void do_cmfd();

        // Vector containing the time that each eigenvalue iteration completed
        // at. Make useful absiccae for convergence plots and the like
        VecF iteration_times_;

    };
}
