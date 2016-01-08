#include "cmfd.hpp"

#include <cmath>
#include <iomanip>
#include <vector>

#include "pugixml.hpp"

#include "mocc-core/error.hpp"
#include "mocc-core/global_config.hpp"

typedef Eigen::Triplet<mocc::real_t> T;
typedef Eigen::SparseMatrix<mocc::real_t> M;

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    CMFD::CMFD( const pugi::xml_node &input, const Mesh *mesh,
            SP_XSMeshHomogenized_t xsmesh):
        mesh_(mesh),
        xsmesh_( xsmesh ),
        n_cell_( mesh->n_pin() ),
        coarse_data_( *mesh, xsmesh->n_group() ),
        is_enabled_( true ),
        fs_( n_cell_ ),
        fs_old_( n_cell_ ),
        source_( n_cell_, xsmesh_.get(), coarse_data_.flux ),
        m_( xsmesh->n_group(), Eigen::SparseMatrix<real_t>(n_cell_, n_cell_) ),
        d_hat_( mesh_->n_surf(), xsmesh_->n_group() ),
        d_tilde_( mesh_->n_surf(), xsmesh_->n_group() ),
        s_hat_( mesh_->n_surf(), xsmesh_->n_group() ),
        s_tilde_( mesh_->n_surf(), xsmesh_->n_group() ),
        k_tol_( 1.0e-4 ),
        psi_tol_( 1.0e-4 ),
        max_iter_( 100 )
    {
        // Set up the structure of the matrix
        std::vector< T > structure;
        for( int i=0; i<n_cell_; i++ ) {
            structure.push_back( T(i, i, 1.0) );
            for( auto d: AllSurfaces ) {
                int n = mesh_->coarse_neighbor( i, d );
                if( n >= 0 ) {
                    // Defining both, though i could do this a little more
                    // efficiently.
                    structure.push_back( T(i, n, 1.0) );
                    structure.push_back( T(n, i, 1.0) );
                }
            }
        }
        for( auto &m: m_ ) {
            m.setFromTriplets( structure.begin(), structure.end() );
            m.makeCompressed();
        }

        // Parse options from the XML, if present
        if( !input.empty() ) {
            // Eigenvalue tolerance
            if( !input.attribute("k_tol").empty() ) {
                k_tol_ = input.attribute("k_tol").as_float(-1.0);
                if( k_tol_ <= 0.0 ) {
                    throw EXCEPT("K tolerance is invalid.");
                }
            }

            // Fission source tolerance
            if( !input.attribute("psi_tol").empty() ) {
                psi_tol_ = input.attribute("psi_tol").as_float(-1.0);
                if( psi_tol_ <= 0.0 ) {
                    throw EXCEPT("Psi tolerance is invalid.");
                }
            }

            // Max iterations
            if( !input.attribute("max_iter").empty() ) {
                max_iter_ = input.attribute("max_iter").as_int(-1);
                if( max_iter_ < 0 ) {
                    throw EXCEPT("Max iterations invalid.");
                }
            }

            // Enabled
            if( !input.attribute("enabled").empty() ) {
                is_enabled_ = input.attribute("enabled").as_bool(true);
            }
        }
        return;
    }


    void CMFD::solve( real_t &k, const ArrayB2 &flux ) {
        // Make sure we have the requisite data
        if( !xsmesh_ ) {
            throw EXCEPT("No XS Mesh data! Need!");
        }

        // Update homogenized cross sections
        xsmesh_->update();
        int ng = xsmesh_->n_group();

        // Set up the linear systems
        this->setup_solve();

        real_t k_old = k;
        real_t tfis = this->total_fission();

        int iter = 0;
        while( true ) {
            iter++;
            // Compute fission source
            fs_old_ = fs_;
            this->fission_source( k );
            real_t tfis_old = tfis;

            for( int group=0; group<ng; group++ ) {
                source_.initialize_group( group );
                source_.fission( fs_, group );
                source_.in_scatter( group );
                source_.scale( mesh_->coarse_volume() );

                this->solve_1g( group );
            }

            tfis = this->total_fission();
            k_old = k;
            k = k * tfis/tfis_old;

            // Convergence check
            real_t psi_err = 0.0;
            for( int i=0; i<(int)fs_.size(); i++ ) {
                real_t e = fs_[i] - fs_old_[i];
                psi_err += e*e;
            }
            psi_err = std::sqrt(psi_err);

            if( ((std::abs(k-k_old) < k_tol_) &&
                (psi_err < psi_tol_)) || (iter > max_iter_) ) {
                break;
            }

        }
        auto flags = cout.flags();
        std::cout << "CMFD : " << iter << " "
                  << std::setprecision(12) << k << " "
                  << std::abs(k-k_old) << std::endl;
        cout.flags(flags);


        // Calculate the resultant currents and store back onto the coarse data
        this->store_currents();

        return;
    }

    void CMFD::solve_1g( int group ) {
        // Not sure exactly how expensive this is. Maybe it would be better to
        // store a collection of BiCGSTAB objects rather than/along with the
        // matrices?
        Eigen::BiCGSTAB< Eigen::SparseMatrix<real_t> > solver;
        solver.compute(m_[group]);
        VectorX x = solver.solve(source_.get());

        // Store the result of the LS solution onto the CoarseData
        ArrayB1 flux_1g = coarse_data_.flux( blitz::Range::all(), group );
        for( int i=0; i<n_cell_; i++ ) {
            flux_1g(i) = x[i];
        }

        return;
    }

    void CMFD::fission_source( real_t k ) {
        int ng = xsmesh_->n_group();

        real_t r_keff = 1.0/k;
        fs_ = 0.0;
        for( int i=0; i<n_cell_; i++ ) {
            for( int ig=0; ig<ng; ig++ ) {
                real_t xsnf = (*xsmesh_)[i].xsmacnf()[ig];
                fs_[i] += r_keff*xsnf*coarse_data_.flux(i, ig);
            }
        }
        return;
    }

    real_t CMFD::total_fission() {
        int ng = xsmesh_->n_group();

        real_t f = 0.0;

        for( int i=0; i<n_cell_; i++ ) {
            for( int ig=0; ig<ng; ig++ ) {
                real_t xsnf = (*xsmesh_)[i].xsmacnf()[ig];
                f += xsnf*coarse_data_.flux(i, ig);
            }
        }

        return f;
    }

    void CMFD::setup_solve() {
        const Mesh::BCArray_t bc = mesh_->boundary_array();
        // Construct the system matrix
        size_t n_surf = mesh_->n_surf();
        int group = 0;
        for( auto &m: m_ ) {
            // Diffusion coefficients
            VecF d_coeff( n_cell_ );
            for( int i=0; i<n_cell_; i++ ) {
                d_coeff[i] = 1.0/(3.0*(*xsmesh_)[i].xsmactr()[group]);
            }

            // Surface diffusivity (d_tilde) and non-linear correction
            // coefficient (d_hat)
            // There are lots of options to optimize this, mostly algebraic
            // simplifications, but this is very conformal to the canonical
            // formulations of CMFD found in the literature. If this starts
            // taking too much time, optimize.
            ArrayB1 d_tilde = d_tilde_( blitz::Range::all(), group );
            ArrayB1 d_hat = d_hat_( blitz::Range::all(), group );
            ArrayB1 s_tilde = s_tilde_( blitz::Range::all(), group );
            ArrayB1 s_hat = s_hat_( blitz::Range::all(), group );
            for( int is=0; is<(int)n_surf; is++ ) {
                auto cells = mesh_->coarse_neigh_cells(is);
                Normal norm = mesh_->surface_normal(is);

                real_t diffusivity_1;
                if( cells.first > -1 ) {
                    // Real neighbor
                    diffusivity_1 = d_coeff[cells.first] /
                        mesh_->cell_thickness(cells.first, norm);
                } else {
                    // Boundary surface
                    switch( bc[(int)norm][0] ) {
                        case Boundary::REFLECT:
                            diffusivity_1 = 0.0;
                            break;
                        case Boundary::VACUUM:
                            diffusivity_1 = 0.5;
                            break;
                        default:
                            cout << "Boundary: " << bc[(int)norm][0] << endl;
                            throw EXCEPT("Unsupported boundary type.");
                    }
                }

                real_t diffusivity_2;
                if( cells.second > -1 ) {
                    // Real neighbor
                    diffusivity_2 = d_coeff[cells.second] /
                        mesh_->cell_thickness(cells.second, norm);
                } else {
                    // Boundary surface
                    switch( bc[(int)norm][1] ) {
                        case Boundary::REFLECT:
                            diffusivity_2 = 0.0;
                            break;
                        case Boundary::VACUUM:
                            diffusivity_2 = 0.5;
                            break;
                        default:
                            cout << "Boundary: " << bc[(int)norm][1] << endl;
                            throw EXCEPT("Unsupported boundary type.");
                    }
                }

                d_tilde(is) = 2.0*diffusivity_1*diffusivity_2 /
                    (diffusivity_1 + diffusivity_2);

                bool have_data = norm == Normal::Z_NORM ?
                    coarse_data_.has_axial_data() :
                    coarse_data_.has_radial_data();
                if( have_data ) {
                    real_t j = coarse_data_.current( is, group );
                    real_t flux_l = cells.first >= 0 ?
                        coarse_data_.flux(cells.first, group) :
                        0.0;
                    real_t flux_r = cells.second >= 0 ?
                        coarse_data_.flux(cells.second, group) :
                        0.0;
                    d_hat(is) = ( j + d_tilde(is)*(flux_r - flux_l)) /
                        (flux_l + flux_r);
                } else {
                    d_hat(is) = 0.0;
                }

            }

            // put values into the matrix. Optimal access patterns in sparse
            // matrix representations are not obvious, so the best way is to
            // iterate through the matrix linearly and act according to the
            // indices (i.e.  row/col) that we get for each coefficient.
            for( int k=0; k<m.outerSize(); k++ ) {
                for( M::InnerIterator it(m, k); it; ++it ) {
                    auto i = it.row();
                    auto j = it.col();
                    if( i == j ) {
                        // Diagonal element
                        real_t xsrm = (*xsmesh_)[i].xsmacrm()[group];
                        real_t v = mesh_->coarse_volume( i ) * xsrm;
                        for( auto is: AllSurfaces ) {
                            size_t surf = mesh_->coarse_surf( i, is );
                            real_t a = mesh_->coarse_area( i, is );
                            // Switch sign of D-hat if necessary
                            real_t d_hat_ij = d_hat(surf);
                            if( (is == Surface::WEST) ||
                                (is == Surface::SOUTH) ||
                                (is == Surface::BOTTOM) )
                            {
                                d_hat_ij = -d_hat_ij;
                            }

                            v += a*(d_tilde(surf) + d_hat_ij);
                        }
                        it.valueRef() = v;
                    } else {
                        // off-diagonal element
                        auto pair = mesh_->coarse_interface(i, j);
                        real_t a = mesh_->coarse_area( i, pair.second );
                        size_t surf = pair.first;
                        real_t d_hat_ij = d_hat(surf);
                        // Switch sign of D-hat if necessary
                        if( (pair.second == Surface::WEST) ||
                            (pair.second == Surface::SOUTH) ||
                            (pair.second == Surface::BOTTOM) )
                        {
                            d_hat_ij = -d_hat_ij;
                        }

                        it.valueRef() = a*(d_hat_ij - d_tilde(surf));
                    }
                }
            } // matrix element loop
            group++;
        } // group loop
        return;
    }

    void CMFD::store_currents() {
        int n_group = xsmesh_->n_group();
        int n_surf = mesh_->n_surf();
        for( int ig=0; ig<n_group; ig++ ) {
            for( int is=0; is<n_surf; is++ ) {
                auto cells = mesh_->coarse_neigh_cells(is);
                real_t flux_r = cells.second >= 0 ?
                    coarse_data_.flux(cells.second, ig) : 0.0;
                real_t flux_l = cells.first >= 0 ?
                    coarse_data_.flux(cells.first, ig) : 0.0;
                real_t d_hat = d_hat_(is, ig);
                real_t d_tilde = d_tilde_(is, ig);

                real_t current = -d_tilde*(flux_r - flux_l) +
                                   d_hat*(flux_r + flux_l);
                coarse_data_.current(is, ig) = current;
            }
        }
        return;
    }
}
