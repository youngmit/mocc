/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <cmath>

#include "core/angle.hpp"
#include "core/angular_quadrature.hpp"
#include "core/global_config.hpp"

namespace mocc { namespace sn {
    /**
     * This is an abstract class, used to define the interfaces that should be
     * provided by all \ref CellWorker derivatives in an Sn sweep. The whole
     * \ref CellWorker concept hinges on the implementation of the \ref
     * CellWorker::evaluate() function being fast and inline, and therefore,
     * there shouldn't be any instances of this class (i.e. \ref CellWorker
     * pointers) floating around, as the polymorphism will kill performance.
     * Rather, this class is provided so that the compiler can check all of the
     * classes that SHOULD implement the right interfaces before we run into
     * cryptic templating errors.
     */
    class CellWorker {
    public:
        CellWorker( const Mesh &mesh, const AngularQuadrature &ang_quad ):
            mesh_( mesh ),
            plane_size_( mesh.nx()*mesh.ny() )
        {
            return;
        }

        inline void set_group( size_t group ) {
            return;
        }

        /**
         * \brief Configure the \ref CellWorker to sweep cells on the given z
         * position.
         */
        inline void set_z( size_t iz ) {
            tz_ = oz_/mesh_.dz(iz);
        }

        /**
         * \brief Configure the \ref CellWorker to sweep cells on the given y
         * position.
         */
        inline void set_y( size_t iy ) {
            ty_ = oy_/mesh_.dy(iy);
        }

        /**
         * \brief Configure the \ref CellWorker to sweep the given angle.
         */
        virtual inline void set_angle( size_t iang, Angle angle ) {
            iang_ = iang;
            angle_ = angle;
            ox_ = std::abs(angle.ox);
            oy_ = std::abs(angle.oy);
            oz_ = std::abs(angle.oz);
        }

        /**
         * \brief Propagate flux through a single mesh element. Return the
         * node-average flux. 3-D.
         *
         * \param[inout] flux_x the upwind flux in the x-normal direction. At
         * the end, it will be updated to the downwind flux.
         * \param[inout] flux_y the upwind flux in the y-normal direction. At
         * the end, it will be updated to the downwind flux.
         * \param[inout] flux_z the upwind flux in the z-normal direction. At
         * the end, it will be updated to the downwind flux.
         * \param[in] q the node-average source.
         * \param[in] xstr the node-average transport cross section.
         * \param[in] i the index of the cell to treat.
         */
        virtual real_t evaluate( real_t &flux_x, real_t &flux_y, real_t &flux_z,
                               real_t q, real_t xstr, size_t i ) = 0;

        /**
         * \brief Propagate flux through a single mesh element in 2-D
         *
         * \param[inout] flux_x the upwind flux in the x-normal direction. At
         * the end, it will be updated to the downwind flux.
         * \param[inout] flux_y the upwind flux in the y-normal direction. At
         * the end, it will be updated to the downwind flux.
         * \param[in] q the node-average source.
         * \param[in] xstr the node-average transport cross section.
         * \param[in] i the index of the cell to treat.
         */
        virtual real_t evaluate_2d( real_t &flux_x, real_t &flux_y,
                               real_t q, real_t xstr, size_t i ) = 0;


    protected:
        const Mesh &mesh_;
        size_t plane_size_;

        real_t ty_;
        real_t tz_;
        size_t iang_;
        Angle angle_;
        real_t ox_;
        real_t oy_;
        real_t oz_;
    };
} }
