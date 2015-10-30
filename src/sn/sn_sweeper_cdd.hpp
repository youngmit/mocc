#pragma once

#include <memory>

#include "pugixml.hpp"

#include "mocc-core/core_mesh.hpp"

#include "sn/cell_worker.hpp"
#include "sn/correction_data.hpp"
#include "sn/sn_sweeper.hpp"

namespace mocc {
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void set_corrections( CorrectionData *data ) {
            corrections_ = data;
            cell_worker_.set_corrections( data );
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
            
            void set_corrections( CorrectionData *data ) {
                corrections_ = data;
            }

            inline real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z,
                                 real_t q, real_t xstr, size_t i )
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

                real_t psi = q + 2.0*(tx*flux_x + ty_*flux_y + tz_*flux_z );
                psi /= tx/gx + ty_/gy + 2.0*tz_ + xstr;
//std::cout << flux_x  << " " << flux_y << " " << flux_z << std::endl;
//std::cout << tx << " " << ty_ << " " << tz_ << std::endl;
//std::cout << ax << " " << ay << " " << b << " " << q << " " << xstr << " " << psi << std::endl;


                flux_x = (psi - gx*flux_x) / gx;
	    		flux_y = (psi - gy*flux_y) / gy;
	    		flux_z = 2.0*psi - flux_z;

                return psi;
            }
        private:
            const AngularQuadrature &ang_quad_;

            const CorrectionData* corrections_;

            size_t iang_alpha_;

            size_t group_;

        };

        CellWorker_CDD cell_worker_;
        const CorrectionData *corrections_;
        std::unique_ptr<const CorrectionData> my_corrections_;

        template <typename CurrentWorker>
        void sweep_cdd( int group );
    };
}
