#pragma once

#include <memory>

#include "pugixml.hpp"

#include "mocc-core/core_mesh.hpp"
#include "mocc-core/files.hpp"

#include "sn/cell_worker.hpp"
#include "sn/correction_data.hpp"
#include "sn/sn_sweeper.hpp"

namespace mocc { namespace sn {
    /**
     * An extension of \ref sn::CellWorker to propagate flux through an
     * orthogonal mesh region with the corrected diamond difference (CDD)
     * scheme. This class is still virtual, as the
     * sn::CellWorker::evaluate() method can be tailored for different axial
     * treatments.
     */
    class CellWorker_CDD: public sn::CellWorker {
    public:
        CellWorker_CDD( const Mesh &mesh,
                const AngularQuadrature &ang_quad ):
            CellWorker( mesh ),
            ang_quad_( ang_quad )
        {
            return;
        }

        void set_group( size_t group ) {
            group_ = group;
        }

        inline void set_angle( size_t iang, Angle angle ) {
            sn::CellWorker::set_angle( iang, angle );
            iang_alpha_ = iang % (ang_quad_.ndir() / 2);
        }

        void set_corrections( const CorrectionData *data ) {
            corrections_ = data;
        }

    protected:
        const AngularQuadrature &ang_quad_;

        const CorrectionData* corrections_;

        size_t iang_alpha_;

        size_t group_;
    };

    /**
     * An extension of \ref CellWorker_CDD to propagate flux through an
     * orthogonal mesh region with the corrected diamond difference (CDD)
     * in X and Y, with diamond difference in Z.
     */
    class CellWorker_CDD_DD: public CellWorker_CDD {
    public:
        CellWorker_CDD_DD( const Mesh &mesh,
                const AngularQuadrature &ang_quad ):
            CellWorker_CDD( mesh, ang_quad )
        {
            return;
        }

        inline real_t evaluate( real_t &flux_x, real_t &flux_y,
                real_t &flux_z, real_t q, real_t xstr, size_t i )
        {
            size_t ix = i % mesh_.nx();
            real_t tx = ox_/mesh_.dx(ix);

            real_t ax = corrections_->alpha( i, iang_alpha_, group_,
                    Normal::X_NORM);
            real_t ay = corrections_->alpha( i, iang_alpha_, group_,
                    Normal::Y_NORM);
            real_t b = corrections_->beta( i, iang_alpha_, group_ );

            real_t gx = ax*b;
            real_t gy = ay*b;

            real_t psi = q + 2.0*(tx * flux_x +
                                  ty_* flux_y +
                                  tz_* flux_z );
            psi /= tx/gx + ty_/gy + 2.0*tz_ + xstr;

            flux_x = (psi - gx*flux_x) / gx;
            flux_y = (psi - gy*flux_y) / gy;
            flux_z = 2.0*psi - flux_z;

            return psi;
        }
    };

    /**
     * An variant of \ref CellWorker_CDD to propagate flux through an
     * orthogonal mesh region with the corrected diamond difference (CDD)
     * scheme in X and Y, with FW difference in Z.
     */
    class CellWorker_CDD_FW: public CellWorker_CDD {
    public:
        CellWorker_CDD_FW( const Mesh &mesh,
                const AngularQuadrature &ang_quad ):
            CellWorker_CDD( mesh, ang_quad )
        {
            return;
        }

        inline real_t evaluate( real_t &flux_x, real_t &flux_y,
                real_t &flux_z, real_t q, real_t xstr, size_t i )
        {
            size_t ix = i % mesh_.nx();
            real_t tx = ox_/mesh_.dx(ix);

            real_t ax = corrections_->alpha( i, iang_alpha_, group_,
                    Normal::X_NORM);
            real_t ay = corrections_->alpha( i, iang_alpha_, group_,
                    Normal::Y_NORM);
            real_t b = corrections_->beta( i, iang_alpha_, group_ );

            real_t gx = ax*b;
            real_t gy = ay*b;

            real_t psi = q + 2.0*(tx * flux_x +
                                  ty_* flux_y ) +
                                  tz_* flux_z ;
            psi /= tx/gx + ty_/gy + tz_ + xstr;

            flux_x = (psi - gx*flux_x) / gx;
            flux_y = (psi - gy*flux_y) / gy;
            flux_z = psi;

            return psi;
        }
    };

    template <class Worker>
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input,
                const CoreMesh &mesh):
            SnSweeper( input, mesh ),
            cell_worker_( mesh, ang_quad_ ),
            corrections_( nullptr )
        {
            LogFile << "Constructing a CDD Sn sweeper" << std::endl;

            if( !input.child("data").empty() ) {
                LogFile << "Located auxiliary data specification." << std::endl;
                if( !input.child("data").attribute("type").empty() ) {
                    std::string data_type =
                        input.child("data").attribute("type").value();
                    if( data_type == "default" ) {
                        LogFile << "Generating default values for correction "
                            "factors." << std::endl;
                        my_corrections_.reset( new CorrectionData( n_reg_,
                                    ang_quad_.ndir(), n_group_) );
                        this->set_corrections( my_corrections_.get() );
                    } else {
                        throw EXCEPT( "Unrecognized data type specified for Sn CDD "
                                "sweeper." );
                    }
                } else {
                    throw EXCEPT( "The <data> tag for an Sn sweeper must have a "
                            "type attrubute." );
                }
            }
        }

        /**
         * \brief Associate the sweeper with a set of correction data.
         */
        void set_corrections( const CorrectionData *data ) {
            // only re-assign the corrections if they are not internally
            // assigned.
            /**
             * \todo It is nice to be able to use default values (0.5) for the
             * corrections, but doing it this way doubles the memory use, since
             * the internally-allocated corrections are stored along with those
             * used by the 2D3D sweeper. There are several ways around this, but
             * i need to decide which way to go.
             */
            if( ( my_corrections_.get() && (data == my_corrections_.get()) ) ||
                !my_corrections_.get() )
            {
                corrections_ = data;
                cell_worker_.set_corrections( data );
            } else {
                LogFile << "CDD sweeper bypassing correction factor assignment "
                    "since they are internally assigned." << std::endl;
            }
        }

        /**
         * \brief Re-assign the angular quadrature.
         */
        void set_ang_quad( AngularQuadrature ang_quad ) {
            ang_quad_ = ang_quad;
            return;
        }

        void sweep( int group ) {
            // Make sure we have correction factors
            if( !corrections_ ) {
                throw EXCEPT( "CDD sweeper doesn't have any correction data. "
                        "Try adding <data type=\"default\"/> in the input file."
                        );
            }
            cell_worker_.set_group( group );

            // Store the transport cross section somewhere useful
            for( auto &xsr: *xs_mesh_ ) {
                real_t xstr = xsr.xsmactr()[group];
                for( auto &ireg: xsr.reg() ) {
                    xstr_[ireg] = xstr;
                }
            }

            // For now this is doing a copy. At some point it might be nice to
            // handle as a reference slice.
            flux_1g_ = flux_(blitz::Range::all(), group);

            // Perform inner iterations
            for( unsigned int inner=0; inner<n_inner_; inner++ ) {
                // Set the source (add upscatter and divide by 4PI)
                source_->self_scatter( group, flux_1g_, q_ );

                if( inner == n_inner_-1 && coarse_data_ ) {
                    // Wipe out the existing currents
                    /// \todo get away from having to manually specify template
                    /// parameters here.
                    coarse_data_->current( blitz::Range::all(), group ) = 0.0;
                    this->sweep_1g<sn::Current, CellWorker_CDD>( group,
                            cell_worker_ );
                    coarse_data_->set_has_data( true );
                    coarse_data_->current( blitz::Range::all(), group ) = 0.0;
                } else {
                    this->sweep_1g<sn::NoCurrent, CellWorker_CDD>( group,
                            cell_worker_ );
                }
            }

            flux_( blitz::Range::all(), group ) = flux_1g_;
            return;
        }

    private:
        Worker cell_worker_;
        std::unique_ptr<const CorrectionData> my_corrections_;
        const CorrectionData *corrections_;
    };
} }
