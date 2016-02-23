#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "core/angular_quadrature.hpp"
#include "core/blitz_typedefs.hpp"
#include "core/coarse_data.hpp"
#include "core/eigen_interface.hpp"
#include "core/global_config.hpp"
#include "core/output_interface.hpp"
#include "core/source.hpp"
#include "core/source_factory.hpp"
#include "core/source_isotropic.hpp"
#include "core/xs_mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"

namespace mocc{
    class TransportSweeper: public HasOutput {
    public:
        TransportSweeper( const pugi::xml_node& input, const CoreMesh& mesh ):
            core_mesh_( &mesh ),
            xs_mesh_( new XSMesh(mesh) ),
            n_reg_( mesh.n_reg() ),
            n_group_( xs_mesh_->n_group() ),
            flux_( n_reg_, n_group_ ),
            flux_old_( n_reg_, n_group_ ),
            vol_( n_reg_ ),
            ang_quad_( input.child("ang_quad") ),
            coarse_data_(nullptr)
        {
            return;
        }

        TransportSweeper( const pugi::xml_node &input ):
            ang_quad_( input.child("ang_quad") ),
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
         * \brief Return a vector containing the pin-homogenizes multi-group
         * scalar flux. The values in the vector are ordered group-major.
         */
        ArrayB2 get_pin_flux() const;

        /**
         * \brief Return a const reference to the Sweeper's \ref
         * AngularQuadrature
         */
        const AngularQuadrature& ang_quad() const {
            return ang_quad_;
        }

        /**
         * Produce pin-homogenized scalar flux for the specified group and store
         * in the passed vector.
         */
        virtual void get_pin_flux_1g( int ig, ArrayB1& flux ) const = 0;

        /**
         * \brief Project a single-group pin mesh-homogenized flux to the fine
         * mesh. Return the residual.
         */
        virtual real_t set_pin_flux_1g( int group, const ArrayB1 &pin_flux )
            = 0;

        /**
         * \brief Project a multi-group pin mesh-homogenized flux to the fine
         * mesh. Return the residual.
         */
        real_t set_pin_flux( const ArrayB2 &pin_flux ) {
            real_t e = 0.0;
            for( int ig=0; ig<(int)n_group_; ig++) {
                ArrayB1 flux_1g = pin_flux(blitz::Range::all(), ig);
                real_t e_g = this->set_pin_flux_1g( ig, flux_1g );
                e += e_g*e_g;
            }
            return std::sqrt( e );
        }

        /**
         * Return a const reference to the MG flux
         */
        const ArrayB2& flux() const {
            return flux_;
        }

        /**
         * Return a reference to the MG flux
         */
        ArrayB2& flux() {
            return flux_;
        }

        /**
         * Given the current estimate of a system eigenvalue, calculate the
         * group-independent fission source and store in the passed array
         */
        virtual void calc_fission_source( real_t k,
                ArrayB1& fission_source) const;

        /**
         * \brief Construct and return a source object which conforms to the
         * sweeper.
         *
         * For now, default to the isotropic MoC Source type, \ref
         * SourceIsotropic.
         */
        virtual UP_Source_t create_source( const pugi::xml_node &input ) const {
            try {
                auto source = SourceFactory( input, n_reg_, xs_mesh_.get(),
                    this->flux() );

                return source;
            }
            catch( Exception e ) {
                std::cerr << e.what() << std::endl;
                throw EXCEPT("Failed to create source.");
            }
        }

        /**
         * Return a shared pointer to a homogenized XS mesh. This is
         * polymorphic, because some sweepers already operate on a homogenized
         * mesh, and there is no need to generate a new one.
         */
        virtual SP_XSMeshHomogenized_t get_homogenized_xsmesh() = 0;

        int n_reg() const {
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
            assert( ig < (int)n_group_ );
            assert( ireg < (int)n_reg_ );

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

        /**
         * Associate the sweeper with a source. This is usually done by
         * something like the FixedSourceSolver.
         */
        virtual void assign_source( Source* source) {
            assert( source != nullptr );
            source_ = source;
        }

        /**
         * \brief Store the current flux as the old flux
         */
        virtual void store_old_flux() {
            flux_old_ = flux_;
            return;
        }

        /**
         * \brief Compute a flux residual between the current state of the flux
         * and the old flux. Defaults to L-2 norm.
         */
        virtual real_t flux_residual() const;

        /**
         * \brief Compute the total fission source based on the current state of
         * the flux
         */
        virtual real_t total_fission( bool old=false ) const;

        /**
         * \brief Return a 3-D array containing normalized pin powers
         *
         * The nature of the normalization is somewhat up in the air. There are
         * different ways to do this. For instance, in the case of non-uniform
         * plane thicknesses and uniform power distribution, should a thicker
         * plane have a greater normalized pin power than a shorter plane? If
         * the pin volumes are not uniform (e.g. annular fuel), should a smaller
         * pin have less normalized power. Generally we would say "no" to the
         * former and "yes" to the latter, but this is pretty arbitrary, so...
         *
         * \note for now the normalization used is really simple! All values are
         * normalized uniformly such that the sum of all powers equals the
         * number of elements in the array.
         */
        virtual ArrayB3 pin_powers() const;

        /**
         * \brief Return a const reference to the region volumes
         */
        const ArrayF &volumes() const {
            if( (int)vol_.size() != n_reg_ ) {
                throw EXCEPT("Volume array dimensions are wrong.");
            }
            return vol_;
        }

    protected:
        const CoreMesh *core_mesh_;

        SP_XSMesh_t xs_mesh_;

        int n_reg_;
        int n_group_;

        Source* source_;

        // Multi-group scalar flux
        ArrayB2 flux_;

        // Previous value of the MG scalar flux
        ArrayB2 flux_old_;

        // Region volumes. In a 3-D sweeper this is the true volume, while in a
        // 2-D sweeper, this is actually surface area.
        ArrayF vol_;

        AngularQuadrature ang_quad_;

        // Reference to the CoarseData object that should be used to store
        // coarse mesh values. This is passed in from above.
        CoarseData *coarse_data_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
