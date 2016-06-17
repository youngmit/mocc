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

#include "eigen_solver.hpp"

#include <iomanip>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"

const static int out_w = 14;

namespace mocc{
    EigenSolver::EigenSolver( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        fss_( input, mesh ),
        fission_source_( fss_.n_reg() ),
        fission_source_prev_( fss_.n_reg() ),
        min_iterations_( 0 )
    {
        LogFile << "Initializing Eigenvalue solver..." << std::endl;


        if( input.empty() ) {
            throw EXCEPT("No input specified for the eigenvalue solver.");
        }

        // grab the convergence constraints from the XML
        int in_int = 0;
        real_t in_float = 0.0;

        // K tolerance
        in_float = input.attribute("k_tol").as_float(-1.0);
        if( in_float <= 0.0 ) {
            throw EXCEPT("Invalid k tolerance.");
        }
        tolerance_k_ = in_float;

        // Psi tolerance
        in_float = input.attribute("psi_tol").as_float(-1.0);
        if( in_float <= 0.0 ) {
            throw EXCEPT("Invalid psi tolerance.");
        }
        tolerance_psi_ = in_float;

        // Max iterations
        in_int = input.attribute("max_iter").as_int(-1);
        if( in_int < 0 ) {
            throw EXCEPT("Invalid number of maximum iterations.");
        }
        max_iterations_ = in_int;

        // Min iterations
        if( !input.attribute("min_iter").empty() ) {
            in_int = input.attribute("min_iter").as_int(-1);
            if( (in_int < 0) || (in_int > (int)max_iterations_) ) {
                throw EXCEPT("Invalid number of minimum iterations.");
            }
            min_iterations_ = in_int;
        }

        // Count the number of fissile mesh regions
        n_fissile_regions_ = 0;
        for( const auto &xsr: fss_.sweeper()->xs_mesh() ) {
            bool has_fission = false;
            for( int ig=0; ig<xsr.n_group(); ig++ ) {
                if( xsr.xsmacnf(ig) > 0.0 ) {
                    has_fission = true;
                    break;
                }
            }
            if( has_fission ) {
                n_fissile_regions_ += xsr.reg().size();
            }
        }

        // CMFD acceleration
        bool do_cmfd = input.attribute("cmfd").as_bool(false);
        if( do_cmfd ) {
            // construct the CMFD solver using the mesh from the transport
            // sweeper
            cmfd_.reset( new CMFD( input.child("cmfd"), (Mesh*)&mesh,
                        fss_.sweeper()->get_homogenized_xsmesh() ) );
            // Associate the sweeper with the coarse data from the CMFD solver
            CoarseData * const cd = cmfd_->get_data();
            fss_.sweeper()->set_coarse_data( cd );
        }

        LogFile << "Done initializing Eigenvalue solver." << std::endl;

        return;
    }

    // Perform a full-blown eigenvalue solve. Start with a guess for the
    // fission source (flat), start doing power iteration. Once we have that
    // working, we will start factoring in CMFD and other fancy tricks
    void EigenSolver::solve() {
        keff_ = 1.0;
        keff_prev_ = 1.0;

        // initialize the fixed source solver and calculation the initial
        // fission source
        fss_.initialize();

        // Hand a reference to the fission source to the fixed source solver
        fss_.set_fission_source(&fission_source_);

        error_k_ = tolerance_k_; // K residual
        error_psi_ = tolerance_psi_; // L-2 norm of the fission source residual

        fss_.sweeper()->calc_fission_source(keff_, fission_source_);

        LogScreen << std::setw(out_w) << "Time"
                  << std::setw(out_w) << "Iter."
                  << std::setw(out_w) << "k"
                  << std::setw(out_w) << "k error"
                  << std::setw(out_w) << "psi error" << std::endl;

        for( size_t n_iterations=0; n_iterations < max_iterations_;
                n_iterations++ )
        {
            this->step();

            // Check for convergence
            error_k_ = fabs(keff_-keff_prev_);

            // use the old fission source to store the difference between new
            // and old, since we will be filling it with new in the next
            // iteration anyways.
            const auto & vol = fss_.sweeper()->volumes();
            real_t efis = 0.0;
            for( int i=0; i<(int)fission_source_.size(); i++ ) {
                real_t e = (fission_source_(i)-fission_source_prev_(i))*vol[i];
                efis += e*e;
            }
            error_psi_ = std::sqrt( efis/n_fissile_regions_ );

            convergence_.push_back(
                    ConvergenceCriteria(keff_, error_k_, error_psi_) );

            iteration_times_.push_back(RootTimer.time());

            this->print( n_iterations+1, convergence_.back() );

            if( n_iterations >= max_iterations_ ) {
                LogScreen << "Maximum number of iterations reached!"
                          << std::endl;
                break;
            }

            if( (error_k_ < tolerance_k_) && (error_psi_ < tolerance_psi_ ) &&
                (n_iterations >= min_iterations_) ) {
                LogScreen << "Convergence criteria satisfied!" << std::endl;
                break;
            }
        }
    } // solve()

    void EigenSolver::step() {

        if( cmfd_ && cmfd_->is_enabled() ) {
            this->do_cmfd();
        }

        // Perform a group sweep with the FSS
        fss_.sweeper()->calc_fission_source(keff_, fission_source_);
        // Store the old fission source
        fission_source_prev_ = fission_source_;
        fss_.step();

        // Get the total fission sources
        real_t tfis1 = fss_.sweeper()->total_fission(false);
        real_t tfis2 = fss_.sweeper()->total_fission(true);

        // update estimate for k
        keff_prev_ = keff_;
        keff_ = keff_ * tfis1/tfis2;

        // update the fission source
        fss_.sweeper()->calc_fission_source(keff_, fission_source_);

        return;
    }

    void EigenSolver::print( int iter, ConvergenceCriteria conv ) {
        LogScreen << std::setw(out_w) << std::fixed << std::setprecision(5)
             << RootTimer.time() << std::setw(out_w)
             << iter << conv << std::endl;
        return;
    }

    void EigenSolver::do_cmfd() {
        assert(cmfd_);
        assert(cmfd_->is_enabled());
        // push homogenized flux onto the coarse mesh, solve, and pull it
        // back.
        cmfd_->coarse_data().flux = fss_.sweeper()->get_pin_flux();

        // Set the convergence criteria for this solve, there are a few ways we
        // can do this

        CMFDConvergence conv = CMFDConvergence::FIXED;
        switch( conv ) {
        case CMFDConvergence::FIXED:
            cmfd_->set_k_tolerance(tolerance_k_/100.0);
            cmfd_->set_psi_tolerance(tolerance_psi_/100.0);
            break;
        case CMFDConvergence::FLOAT:
            real_t k_tol = std::max(error_k_/1000.0, tolerance_k_/10.0);
            cmfd_->set_k_tolerance( k_tol );

            real_t psi_tol = std::max(error_psi_/1000.0, tolerance_psi_/10.0);
            cmfd_->set_psi_tolerance( psi_tol );
            break;
        }
        cmfd_->solve( keff_ );
        fss_.sweeper()->set_pin_flux( cmfd_->flux() );
        fss_.sweeper()->update_incoming_flux();
        return;
    }

    void EigenSolver::output( H5Node &file ) const {
        VecF k;
        VecF error_k;
        VecF error_psi;

        for( auto &c: convergence_ ) {
            k.push_back(c.k);
            error_k.push_back(c.error_k);
            error_psi.push_back(c.error_psi);
        }

        VecI dims(1, convergence_.size());

        {
            auto g = file.create_group("convergence");

            g.write( "k", k, dims );
            g.write( "error_k", error_k, dims );
            g.write( "error_psi", error_psi, dims );
            g.write( "iteration_time", iteration_times_ );
            g.write( "abscissae", iteration_times_ );
        }

        fss_.output( file );
        return;
    }

    std::ostream& operator<<( std::ostream &os, ConvergenceCriteria conv) {
        std::ios::fmtflags flags = os.flags();

        os   << std::setw(out_w) << std::fixed
                 << std::setprecision(10) << conv.k
             << std::setw(out_w) << std::scientific
                 << std::setprecision(6) << conv.error_k
             << std::setw(out_w) << std::setiosflags(std::ios::scientific)
                 << std::setprecision(6) << conv.error_psi;

        os.flags(flags);
        return os;
    }
};
