#pragma once

#include <memory>
#include <vector>

#include "coarse_data.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "output_interface.hpp"
#include "source.hpp"
#include "xs_mesh.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc{
    class TransportSweeper: public HasOutput {
    public:
        TransportSweeper():
            coarse_data_(nullptr)
        {
            return;
        }

        TransportSweeper( const CoreMesh& mesh ):
            xs_mesh_( new XSMesh(mesh) ),
            n_reg_( mesh.n_reg() ),
            ng_( xs_mesh_->n_group() ),
            flux_( n_reg_, ng_ ),
            flux_old_( n_reg_, ng_ ),
            vol_( n_reg_, 1 ),
            coarse_data_(nullptr)
        {
            return;
        }

        virtual ~TransportSweeper(){ }

        /**
         * Perform a transport sweep of the passed group.
         */
        virtual void sweep(int group) = 0;

        /** 
         * Initialize the solution variables (scalar, boundary flux, etc.) to
         * reasonable initial guesses.
         */
        virtual void initialize() = 0;

        /**
         * Produce pin-homogenized scalar flux for the specified group and store
         * in the passed vector.
         */
        virtual void get_pin_flux( int ig, VecF& flux ) const = 0;

        /**
         * Given the current estimate of a system eigenvalue, calculate the
         * group-independent fission source and store in the passed array
         */
        virtual void calc_fission_source( float_t k, 
                ArrayX& fission_source) const;

        /**
         * Construct and return a source object which conforms to the sweeper.
         * For now, default to the MoC Source type
         */
        virtual UP_Source_t create_source() const {
            UP_Source_t source( new Source( n_reg_, xs_mesh_.get(), 
                        this->cflux()) );
            return source;
        }

        /**
         * Return a shared pointer to a homogenized XS mesh. This is
         * polymorphic, because some sweepers already operate on a homogenized
         * mesh, and there is no need to generate a new one.
         */
        virtual SP_XSMeshHomogenized_t get_homogenized_xsmesh() = 0;
        
        unsigned int n_reg() const {
            return n_reg_;
        }

        /// Return a reference to the sweeper's XSMesh
        const XSMesh& xs_mesh() const {
            return *(xs_mesh_.get());
        }

        /// Return a shared pointer to the sweeper's XSMesh. Use with caution
        SP_XSMesh_t get_xs_mesh() {
            return xs_mesh_;
        }

        /// Return a reference to the MG flux
        const ArrayX& flux() const {
            return flux_;
        }

        /// Subscript and return a specific flux value
        const float_t flux( unsigned int ig, unsigned int ireg ) const {
            return flux_( ireg, ig );
        }

        /// Return the number of energy groups
        unsigned int n_grp() const {
            return ng_;
        }

        /// Assign a CoarseData object to the sweeper, allowing it to store
        /// currents and such.
        virtual void set_coarse_data( CoarseData *cd ) {
            coarse_data_ = cd;
        }

        /// Return a const reference to the MG flux. This is the same as the
        /// above for now, since im not sure if i want to expose a non-const
        /// reference. Probably will at some point, we will see. It'll be less
        /// refactoring if I start with an explicit const version and use it
        /// wherever I know I won't need mutability.
        const ArrayX& cflux() const {
            return flux_;
        }

        /// Associate the sweeper with a source. This is usually done by
        /// something like the FixedSourceSolver.
        virtual void assign_source( const Source* source) {
            assert( source != nullptr );
            source_ = source;
        }

        /// Store the current flux as the old flux
        virtual void store_old_flux() {
            flux_old_ = flux_;
            return;
        }

        /// Compute the total fission source based on the current state of the
        /// flux
        virtual float_t total_fission( bool old=false ) const;

        /// Homogenize flux and group constants to a CoarseData object
        virtual void homogenize( CoarseData &data ) const = 0;

    protected:
        SP_XSMesh_t xs_mesh_;

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

        // Reference to the CoarseData object that should be used to store
        // coarse mesh values. This is passed in from above.
        CoarseData *coarse_data_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
