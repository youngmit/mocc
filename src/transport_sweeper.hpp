#pragma once

#include <memory>
#include <vector>

#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "xs_mesh.hpp"
#include "source.hpp"

namespace mocc{
    class TransportSweeper{
    public:
        TransportSweeper( const CoreMesh& mesh ):
            mesh_( mesh ),
            xs_mesh_( mesh ),
            n_reg_( mesh.n_reg() ),
            ng_( xs_mesh_.n_grp() ),
            flux_( ng_, n_reg_ )
        {
            return;
        }
        virtual ~TransportSweeper(){ }
        virtual void sweep(int group) = 0;
        virtual void initialize() = 0;

        // Given the current estimate of a system eigenvalue, calculate the
        // group-independent fission source and store in the passed array
        virtual void calc_fission_source( float_t k, 
                MatrixX& fission_source) const = 0;
        
        unsigned int n_reg() const {
            return n_reg_;
        }

        // Return a reference to the sweeper's XSMesh
        const XSMesh& xs_mesh() const {
            return xs_mesh_;
        }

        // Return a reference to the MG flux
        const MatrixX& flux() const {
            return flux_;
        }

        // Subscript and return a specific flux value
        const float_t flux( unsigned int ig, unsigned int ireg ) const {
            return flux_(ig, ireg);
        }

        // Return the number of energy groups
        unsigned int n_grp() const {
            return ng_;
        }

        // Return a const reference to the MG flux. This is the same as the
        // above for now, since im not sure if i want to expose a non-const
        // reference. Probably will at some point, we will see. It'll be less
        // refactoring if I start with an explicit const version and use it
        // wherever I know I won't need mutability.
        const MatrixX& cflux() const {
            return flux_;
        }
    protected:
        const CoreMesh& mesh_;
        XSMesh xs_mesh_;

        unsigned int n_reg_;
        unsigned int ng_;

        const Source* source_;

        // Multi-group scalar flux
        MatrixX flux_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
