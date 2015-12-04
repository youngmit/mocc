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
            CellWorker( mesh, ang_quad ),
            ang_quad_( ang_quad ),
            corrections_( nullptr )
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

        /**
         * \brief Associate the internal reference to correction data.
         *
         * Any existing data will get kicked off. Since this uses
         * std::shared_ptr, if the sweeper has the only reference to any data
         * that gets replaced, we should expect the old data to be destroyed.
         * Usually what we want, but be careful.
         */
        void set_corrections( std::shared_ptr<const CorrectionData> data ) {
            corrections_ = data;
        }

        /** 
         * \copydoc CellWorker::evaluate_2d()
         *
         * Since the variants of the CDD worker are all for different axial
         * treatments, the 2-D version of \ref CellWorker::evaluate() can live
         * here.
         */
        inline real_t evaluate_2d( real_t &flux_x, real_t &flux_y, real_t q, 
                real_t xstr, size_t i )
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
                                  ty_* flux_y );
            psi /= tx/gx + ty_/gy + xstr;

            flux_x = (psi - gx*flux_x) / gx;
            flux_y = (psi - gy*flux_y) / gy;

            return psi;
        }

    protected:
        const AngularQuadrature &ang_quad_;

        std::shared_ptr<const CorrectionData> corrections_;

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
     * A variant of \ref CellWorker_CDD to propagate flux through an orthogonal
     * mesh region with the corrected diamond difference (CDD) scheme in X and
     * Y, with FW difference in Z.
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
} }
