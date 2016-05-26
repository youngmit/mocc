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

#include <vector>

#include "geometry/direction.hpp"
#include "geometry/geom.hpp"
#include "global_config.hpp"
#include "lattice.hpp"

namespace mocc {
    class Plane {
    public:
        Plane(const std::vector<const Lattice*> &lattices, size_t nx,
                size_t ny);

        const Lattice& at(size_t ix, size_t iy) const {
            assert(ix >= 0);
            assert(iy >= 0);
            assert(ix < nx_);
            assert(iy < ny_);
            return *(lattices_[ix + nx_*iy]);
        }

        /**
         *
         * \brief Given a \ref Point2 in core-local coordinates, return a const
         * pointer to the corresponding \ref PinMesh.
         *
         * \param[in,out] p a Point2 in core-local coordinates. Will be modified
         * (see below).
         * \param[in,out] first_reg the first FSR index of the Plane. Will be
         * updated to the first region of the \ref Pin in which the pin resides.
         * \param[in] dir optional \ref Direction to use to disambiguate when
         * the \p p lies directly on a border.
         */
        const PinMesh* get_pinmesh( Point2 &p, int &first_reg,
                Direction dir=Direction()) const;

        /**
         * \brief Return a const pointer to the \ref PinMesh that is at the
         * passed \ref Position.
         */
        const PinMesh* get_pinmesh( Position pos ) const;

        /**
         * Return the number of solution mesh regions in the \ref Plane
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * Return the number of XS Mesh regions in the Plane
         */
        size_t n_xsreg() const {
            return n_xsreg_;
        }

        /**
         * \brief Return a vector containing the FSR volumes
         */
        VecF vols() const {
            VecF vols;
            for( auto &lat: lattices_ ) {
                for( auto &pin: *lat ) {
                    vols.insert(vols.end(), pin->vols().begin(),
                            pin->vols().end());
                }
            }

            return vols;
        }

        /**
         * \brief Return the position of a pin, given its index.
         */
        Position pin_position( size_t ipin ) const;

        /**
         * \brief Return the number of \ref Pin s in this \ref Plane marked as
         * fuel.
         */
        int n_fuel() const {
            return n_fuel_;
        }

    private:
        /**
         * Number of lattices in the x direction
         */
        unsigned nx_;
        /**
         * Number of lattices in the y direction
         */
        unsigned ny_;

        size_t n_reg_;
        size_t n_xsreg_;

        /**
         * Locations of \ref Lattice interfaces
         */
        VecF hx_;
        VecF hy_;

        /**
         * Local list of \ref Lattice pointers
         */
        std::vector<const Lattice*> lattices_;

        /**
         * List of the starting FSR index for each \ref Lattice in the plane
         */
        VecI first_reg_lattice_;

        int n_fuel_;
    };
}
