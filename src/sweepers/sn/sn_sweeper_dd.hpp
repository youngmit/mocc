#pragma once

#include "core/mesh.hpp"

#include "sn/cell_worker.hpp"
#include "sn/sn_sweeper.hpp"

namespace mocc { namespace sn {
    /**
     * A simple class extending \ref sn::CellWorker to perform the algebraic
     * work needed to propagate flux through an orthogonal mesh cell using
     * the diamond difference scheme.
     */
    class CellWorker_DD: public sn::CellWorker {
    public:
        CellWorker_DD( const Mesh &mesh, const AngularQuadrature &ang_quad ):
            CellWorker( mesh, ang_quad )
        {
            return;
        }


        inline real_t evaluate(real_t &flux_x, real_t &flux_y,
                real_t &flux_z, real_t q, real_t xstr, size_t i )
        {
            size_t ix = i % mesh_.nx();
            real_t tx = ox_/mesh_.dx(ix);
            real_t psi = 2.0*( tx * flux_x +
                               ty_* flux_y +
                               tz_* flux_z ) + q;
            psi /= 2.0*(tx + ty_ + tz_) + xstr;

            flux_x = 2.0*psi - flux_x;
            flux_y = 2.0*psi - flux_y;
            flux_z = 2.0*psi - flux_z;

            return psi;
        }

        inline real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q,
                real_t xstr, size_t i )
        {
            size_t ix = i % mesh_.nx();
            real_t tx = ox_/mesh_.dx(ix);
            real_t psi = 2.0*( tx * flux_x +
                               ty_* flux_y ) + q;
            psi /= 2.0*(tx + ty_) + xstr;

            flux_x = 2.0*psi - flux_x;
            flux_y = 2.0*psi - flux_y;

            return psi;
        }
    };
} }
