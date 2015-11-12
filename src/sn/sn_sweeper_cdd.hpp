#pragma once

#include <memory>

#include "pugixml.hpp"

#include "mocc-core/core_mesh.hpp"
#include "mocc-core/files.hpp"

#include "sn/cell_worker.hpp"
#include "sn/correction_data.hpp"
#include "sn/sn_sweeper.hpp"

namespace mocc {
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

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

    private:
        /**
         * An extension of \ref sn::CellWorker to propagate flux through an
         * orthogonal mesh region with the corrected diamond difference (CDD)
         * scheme.
         */
        class CellWorker_CDD: public sn::CellWorker {
        public:
            CellWorker_CDD( const Mesh &mesh, const AngularQuadrature &ang_quad ):
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
                flux_z = psi; // override DD with FW diff

                return psi;
            }
        private:
            const AngularQuadrature &ang_quad_;

            const CorrectionData* corrections_;

            size_t iang_alpha_;

            size_t group_;
        };

        CellWorker_CDD cell_worker_;
        std::unique_ptr<const CorrectionData> my_corrections_;
        const CorrectionData *corrections_;

        template <typename CurrentWorker>
        void sweep_cdd( int group );
    };
}
