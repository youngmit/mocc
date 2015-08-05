#pragma once

#include <memory>
#include <vector>

#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "xs_mesh.hpp"
#include "source.hpp"
#include "output_interface.hpp"

namespace mocc{
    class TransportSweeper: public HasOutput {
    public:
        TransportSweeper( ) { }
        TransportSweeper( const CoreMesh& mesh ):
            mesh_( &mesh ),
            n_reg_( mesh.n_reg() ),
            ng_( xs_mesh_.n_grp() ),
            flux_( n_reg_, ng_ ),
            flux_old_( n_reg_, ng_ ),
            vol_( n_reg_, 1 )
        {
            return;
        }

        virtual ~TransportSweeper(){ }
        virtual void sweep(int group) = 0;
        virtual void initialize() = 0;
        virtual void get_pin_flux( int ig, VecF& flux ) const = 0;

        // Given the current estimate of a system eigenvalue, calculate the
        // group-independent fission source and store in the passed array
        virtual void calc_fission_source( float_t k, 
                ArrayX& fission_source) const = 0;
        
        unsigned int n_reg() const {
            return n_reg_;
        }

        // Return a reference to the sweeper's XSMesh
        const XSMesh& xs_mesh() const {
            return xs_mesh_;
        }

        // Return a reference to the MG flux
        const ArrayX& flux() const {
            return flux_;
        }

        // Subscript and return a specific flux value
        const float_t flux( unsigned int ig, unsigned int ireg ) const {
            return flux_( ireg, ig );
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
        const ArrayX& cflux() const {
            return flux_;
        }

        // Associate the sweeper with a source. This is usually done by
        // something like the FixedSourceSolver.
        void assign_source( Source* source) {
            assert( source != nullptr );
            source_ = source;
        }

        // Return the energy group bounds from the underlying xsmesh
        const VecF& eubounds() const {
            return xs_mesh_.eubounds();
        }

        // Store the current flux as the old flux
        void store_old_flux() {
            flux_old_ = flux_;
            return;
        }

        // Compute the total fission source based on the current state of the
        // flux
        float_t total_fission( bool old=false ) const;
    protected:
        const CoreMesh* mesh_;
        XSMesh xs_mesh_;

        unsigned int n_reg_;
        unsigned int ng_;

        const Source* source_;

        // Multi-group scalar flux
        ArrayX flux_;

        // Previous value of the MG scalar flux
        ArrayX flux_old_;

        // Region volumes. In a 3-D sweeper this is the true volume, while in a
        // 2-D sweeper, this is actually surface area.
        ArrayX vol_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
