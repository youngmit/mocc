#include "eigen_solver.hpp"

#include "error.hpp"

namespace mocc{
    EigenSolver::EigenSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        fss_( input, mesh ),
        fission_source_( fss_.n_reg(), 1 ),
        fission_source_prev_( fss_.n_reg(), 1 )
    {
        if( input.empty() ) {
            Error("No input specified for the eigenvalue solver.");
        }

        // grab the convergence constraints from the XML
        int in_int = 0;
        float_t in_float = 0.0;

        // K tolerance
        in_float = input.attribute("k_tol").as_float(-1.0);
        if( in_float <= 0.0 ) {
            Error("Invalid k tolerance.");
        }
        tolerance_k_ = in_float;

        // Psi tolerance
        in_float = input.attribute("psi_tol").as_float(-1.0);
        if( in_float <= 0.0 ) {
            Error("Invalid psi tolerance.");
        }
        tolerance_psi_ = in_float;

        // Max iterations
        in_int = input.attribute("max_iter").as_int(-1);
        if( in_int < 0 ) {
            Error("Invalid number of maximum iterations.");
        }
        max_iterations_ = in_int;
    }

    // Perform a full-blown eigenvalue solve. Start with a guess for the
    // fission source (flat), start doing power iteration. Once we have that
    // working, we will start factoring in CMFD and other fancy tricks
    void EigenSolver::solve() {
        keff = 1.0;
        keff_prev = 1.0;

        // initialize the fixed source solver and calculation the initial
        // fission source
        fss_.initialize();
        
        // Fill the fission source with zeros. when its copied into the old
        // fission source on the first iteration, we should end up with a pretty
        // large difference, which ought not to satisfy our convergence criteria
        fission_source_.fill(0.0);


        float_t error_k = 1.0; // K residual
        float_t error_psi = 1.0; // L-2 norm of the fission source residual

        unsigned int n_iterations = 0;

        bool done = false;
        while( !done ) {
            
            done = (error_k < tolerance_k_) & (error_psi < tolerance_psi_) & 
                (n_iterations >= max_iterations_ );
        }

    }

    void EigenSolver::step() {

        // Store the old fission source
        fission_source_prev_ = fission_source_;


        fss_.sweeper()->calc_fission_source(keff, fission_source_);
        
    }
    
};
