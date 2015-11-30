#include "eigen_solver.hpp"

#include <iomanip>

#include "error.hpp"

using std::endl;
using std::cout;
using std::cin;

const static int out_w = 14;

namespace mocc{
    EigenSolver::EigenSolver( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        fss_( input, mesh ),
        fission_source_( fss_.n_reg() ),
        fission_source_prev_( fss_.n_reg() )
    {
        if( input.empty() ) {
            Error("No input specified for the eigenvalue solver.");
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

        real_t error_k = 1.0; // K residual
        real_t error_psi = 1.0; // L-2 norm of the fission source residual

        unsigned int n_iterations = 0;

        bool done = false;

        cout << std::setw(out_w) << "Iter."
             << std::setw(out_w) << "k"
             << std::setw(out_w) << "k error"
             << std::setw(out_w) << "psi error" << endl;

        while( !done ) {

            this->step();

            n_iterations++;

            // Check for convergence
            error_k = fabs(keff_-keff_prev_);

            // use the old fission source to store the difference between new
            // and old, since we will be filling it with new in the next
            // iteration anyways.
            error_psi = std::sqrt(
                    std::pow(fission_source_-fission_source_prev_, 2).sum() );

            convergence_.push_back(
                    ConvergenceCriteria(keff_, error_k, error_psi) );

            this->print( n_iterations, convergence_.back() );

            done = ((error_k < tolerance_k_) & (error_psi < tolerance_psi_)) |
                (n_iterations >= max_iterations_ );
        }

    }

    void EigenSolver::step() {
        // Store the old fission source
        fission_source_prev_ = fission_source_;

        if( cmfd_ ) {
            // push homogenized flux onto the coarse mesh
            cmfd_->coarse_data().flux = fss_.sweeper()->get_pin_flux();
            cout << "flux into cmfd:" << endl;
            cout << fss_.sweeper()->get_pin_flux() << endl;;
            cmfd_->solve(keff_, fss_.sweeper()->flux());
            fss_.sweeper()->set_pin_flux( cmfd_->flux() );
        }


        // Perform a group sweep with the FSS
        fss_.sweeper()->calc_fission_source(keff_, fission_source_);
        fss_.step();


        // Get the total fission sources
        real_t tfis1 = fss_.sweeper()->total_fission(false);
        real_t tfis2 = fss_.sweeper()->total_fission(true);

        // update estimate for k
        keff_prev_ = keff_;
        keff_ = keff_ * tfis1/tfis2;
    }

    void EigenSolver::print( int iter, ConvergenceCriteria conv ) {
        cout  << std::setw(out_w) << iter << conv << endl;;
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
