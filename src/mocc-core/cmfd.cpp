#include "cmfd.hpp"

#include <vector>

#include "mocc-core/error.hpp"
#include "mocc-core/global_config.hpp"

typedef Eigen::Triplet<mocc::real_t> T;
typedef Eigen::SparseMatrix<mocc::real_t> M;

namespace mocc {
    CMFD::CMFD( const Mesh *mesh, SP_XSMeshHomogenized_t xsmesh):
        mesh_(mesh),
        xsmesh_( xsmesh ),
        coarse_data_( mesh->n_pin(), mesh->n_surf(), xsmesh->n_group() ),
        b_( mesh->n_pin() ),
        source_( mesh->n_pin(), xsmesh_.get(), coarse_data_.flux )
    {
        // Set up the structure of the matrix
        size_t n_cell = mesh_->n_pin();
        std::vector< T > structure;
        for( size_t i=0; i<n_cell; i++ ) {
            structure.push_back( T(i, i, 1.0) );
            for( auto d: AllSurfaces ) {
                size_t n = mesh_->coarse_neighbor( i, d );
                if( n >= 0 ) {
                    // Defining both, though i could do this a little more
                    // efficiently.
                    structure.push_back( T(i, n, 1.0) );
                    structure.push_back( T(n, i, 1.0) );
                }
            }
        }
        m_.setFromTriplets( structure.begin(), structure.end() );
        m_.makeCompressed();
        return;
    }


    void CMFD::solve( real_t &k, const ArrayB2 &flux ) {
        // Make sure we have the requisite data
        if( !xsmesh_ ) {
            throw EXCEPT("No XS Mesh data! Need!");
        }

        // Update homogenized cross sections
        xsmesh_->update( flux );
        int ng = xsmesh_->n_group();

        real_t k_old = k;
        while( true ) {
            // Compute fission source
            

            for( int group=0; group<ng; group++ ) {
            
            }
        }

        return;
    }

    void CMFD::solve_1g( int group ) {
        const auto bc = mesh_->boundary_array();
        // Construct the system matrix
        size_t n_surf = mesh_->n_surf();
        size_t n_cell = mesh_->n_pin();
        //
        // Diffusion coefficients
        VecF d_coeff( 0.0, mesh_->n_pin() );
        for( size_t i=0; i<mesh_->n_pin(); i++ ) {
            d_coeff[i] = 1.0/(3.0*(*xsmesh_)[i].xsmactr()[group]);
        }

        // Surface diffusivity (d_tilde) and non-linear correction coefficient
        // (d_hat)
        // There are lots of options to optimize this, mostly algebraic
        // simplifications, but this is very conformal to the canonical
        // formulations of CMFD found in the literature. If this starts taking
        // too much time, optimize.
        VecF d_tilde( n_surf, 0.0 );
        VecF d_hat( n_surf, 0.0 );
        for( size_t is=0; is<n_surf; is++ ) {
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
                        diffusivity_1 = 0.5;
                        break;
                    case Boundary::VACUUM:
                        diffusivity_1 = 0.0;
                        break;
                    default:
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
                switch( bc[(int)norm][0] ) {
                    case Boundary::REFLECT:
                        diffusivity_2 = 0.5;
                        break;
                    case Boundary::VACUUM:
                        diffusivity_2 = 0.0;
                        break;
                    default:
                        throw EXCEPT("Unsupported boundary type.");
                }
            }

            d_tilde[is] = 2.0*diffusivity_1*diffusivity_2 / 
                (diffusivity_1 + diffusivity_2);

            real_t j = coarse_data_.current( is, group );
            real_t flux_l = coarse_data_.flux(cells.first, group);
            real_t flux_r = coarse_data_.flux(cells.second, group);
            d_hat[is] = ( j - d_tilde[is]*(flux_r - flux_l)) / 
                (flux_l + flux_r);
        }

        // put values into the matrix. Optimal access patterns in sparce matrix
        // representations are not obvious, so the best way is to iterate
        // through the matrix linearly and act according to the indices (i.e.
        // row/col) that we get for each coefficient.
        for( int k=0; k<m_.outerSize(); k++ ) {
            for( M::InnerIterator it(m_, k); it; ++it ) {
                auto i = it.row();
                auto j = it.col();
                if( i == j ) {
                    // Diagonal element
                    real_t xsrm = (*xsmesh_)[i].xsmacrm()[group];
                    real_t v = mesh_->coarse_volume( i ) * xsrm;
                    for( auto is: AllSurfaces ) {
                        size_t surf = mesh_->coarse_surf( i, is );
                        real_t a = mesh_->coarse_area( i, is );
                        v += a*(d_tilde[surf] + d_hat[surf]);
                    }
                    it.valueRef() = v;
                } else {
                    // off-diagonal element
                    auto pair = mesh_->coarse_interface(i, j);
                    real_t a = mesh_->coarse_area( i, pair.second );
                    size_t surf = pair.first;
                    it.valueRef() = a*(d_hat[surf] - d_tilde[surf]);
                }
            }
        }

        Eigen::BiCGSTAB< Eigen::SparseMatrix<real_t> > solver;
        solver.compute(m_);
        VectorX x = solver.solve(b_);
    }
}
