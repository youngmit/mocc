#pragma once

#include <memory>
#include <vector>

#include "mocc-core/blitz_typedefs.hpp"
#include "mocc-core/coarse_data.hpp"
#include "mocc-core/eigen_interface.hpp"
#include "mocc-core/global_config.hpp"
#include "mocc-core/output_interface.hpp"
#include "mocc-core/source.hpp"
#include "mocc-core/xs_mesh.hpp"
#include "mocc-core/xs_mesh_homogenized.hpp"

namespace mocc{
    class TransportSweeper: public HasOutput {
    public:
        TransportSweeper():
            coarse_data_(nullptr)
        {
            return;
        }

        TransportSweeper( const CoreMesh& mesh ):
            core_mesh_( &mesh ),
            xs_mesh_( new XSMesh(mesh) ),
            n_reg_( mesh.n_reg() ),
            n_group_( xs_mesh_->n_group() ),
            flux_( n_reg_, n_group_ ),
            flux_old_( n_reg_, n_group_ ),
            vol_( n_reg_ ),
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
        virtual void get_pin_flux_1g( int ig, VecF& flux ) const = 0;
        
        /**
         * Return a reference to the MG flux
         */
        const ArrayB2& flux() const {
            return flux_;
        }

        /**
         * Return a vector containing the pin-homogenizes multi-group scalar
         * flux. The values in the vector are ordered group-major.
         */
        VecF get_pin_flux() const;

        /**
         * Given the current estimate of a system eigenvalue, calculate the
         * group-independent fission source and store in the passed array
         */
        virtual void calc_fission_source( real_t k,
                ArrayF& fission_source) const;

        /**
         * Construct and return a source object which conforms to the sweeper.
         * For now, default to the MoC Source type
         */
        virtual UP_Source_t create_source() const {
            UP_Source_t source( new Source( n_reg_, xs_mesh_.get(),
                        this->flux()) );
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

        /**
         * Return a reference to the sweeper's XSMesh
         */
        const XSMesh& xs_mesh() const {
            return *(xs_mesh_.get());
        }

        /**
         * Return a shared pointer to the sweeper's XSMesh. Use with caution
         */
        SP_XSMesh_t get_xs_mesh() {
            return xs_mesh_;
        }

        /**
         * Return a reference to the CoreMesh
         */
        const CoreMesh& mesh() const {
            return *core_mesh_;
        }

        /**
         * Subscript and return a specific flux value
         */
        real_t flux( int ig, int ireg ) const {
            assert( ig < n_group_ );
            assert( ireg < n_reg_ );

            return flux_( ireg, ig );
        }

        /**
         * Return the number of energy groups
         */
        unsigned int n_group() const {
            return n_group_;
        }

        /**
         * Assign a CoarseData object to the sweeper, allowing it to store
         * currents and such.
         */
        virtual void set_coarse_data( CoarseData *cd ) {
            coarse_data_ = cd;
        }

        /// Associate the sweeper with a source. This is usually done by
        /// something like the FixedSourceSolver.
        virtual void assign_source( Source* source) {
            assert( source != nullptr );
            source_ = source;
        }

        /// Store the current flux as the old flux
        virtual void store_old_flux() {
            flux_old_ = flux_;
            return;
        }

        /**
         * Compute a flux residual between the current state of the flux and the
         * old flux. Defaults to L-2 norm.
         */
        virtual real_t flux_residual() const;

        /**
         * Compute the total fission source based on the current state of the
         * flux
         */
        virtual real_t total_fission( bool old=false ) const;

        /**
         * Homogenize flux and group constants to a CoarseData object
         */
        virtual void homogenize( CoarseData &data ) const = 0;

        /**
         * \brief Project a pin-mesh flux to the fine mesh. Return the residual.
         */
        virtual real_t set_pin_flux_1g( int group, const VecF &pin_flux ) = 0;

    protected:
        const CoreMesh *core_mesh_;

        SP_XSMesh_t xs_mesh_;

        unsigned int n_reg_;
        unsigned int n_group_;

        Source* source_;

        // Multi-group scalar flux
        ArrayB2 flux_;

        // Previous value of the MG scalar flux
        ArrayB2 flux_old_;

        // Region volumes. In a 3-D sweeper this is the true volume, while in a
        // 2-D sweeper, this is actually surface area.
        ArrayF vol_;

        // Reference to the CoarseData object that should be used to store
        // coarse mesh values. This is passed in from above.
        CoarseData *coarse_data_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
