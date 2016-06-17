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

#include <iosfwd>
#include <memory>
#include <map>
#include <vector>

#include "core/geometry/direction.hpp"
#include "core/geometry/geom.hpp"
#include "core/global_config.hpp"
#include "core/pin.hpp"
#include "core/pin_mesh_base.hpp"
#include "core/pugifwd.hpp"

namespace mocc {
    class Lattice {
    public:
        Lattice( const pugi::xml_node &input,
                 const std::map<int, UP_Pin_t> &pins );

        size_t id() const {
            return id_;
        }

        /**
         * \brief Number of pins in the x direction
         */
        size_t nx() const {
            return nx_;
        }

        /**
         * \brief Number of pins in the y direction
         */
        size_t ny() const {
            return ny_;
        }

        /**
         * \brief Total number of pins in the lattice
         */
        size_t n_pin() const {
            return pins_.size();
        }

        /**
         * \brief Return the size of the lattice along the x dimension
         */
        real_t hx() const {
            return hx_;
        }

        /**
         * \brief Return the size of the lattice along the y dimension
         */
        real_t hy() const {
            return hy_;
        }

        /**
         * \brief Return a const reference to the \ref Pin located at the given
         * location of the \ref Lattice.
         */
        const Pin& at( size_t x, size_t y ) const {
            assert( (0 <= x) & (x < nx_) );
            assert( (0 <= y) & (y < ny_) );
            return *pins_[y*nx_ + x];
        }

        /**
         * \brief Return an iterator to the first \ref Pin* in the \ref Lattice
         */
        std::vector<Pin*>::const_iterator begin() const {
            return pins_.cbegin();
        }

        /**
         * \brief Return an iterator past the last Pin* in the lattice
         */
        std::vector<Pin*>::const_iterator end() const {
            return pins_.cend();
        }

        /**
         * \brief Return a const reference to the underlying hx_vec_ array
         */
        const VecF& hx_vec() const {
            return hx_vec_;
        }

        /**
         * \brief Return a const reference to the underlying hy_vec_ array
         */
        const VecF& hy_vec() const {
            return hy_vec_;
        }

        /**
         * \brief Return the total number of regions in the lattice
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * \brief Return the total number of XS regions in the lattice
         */
        size_t n_xsreg() const {
            return n_xsreg_;
        }

        /**
         * \brief Return a const pointer to the \ref PinMesh object located at
         * the provided point, incrementing the passed in first_reg by the pin's
         * first-region offset.
         *
         * \param[in,out] p a \ref Point2 in lattice-local coordinates to look
         * up. Will be updated to the location of the pin origin.
         * \param[in,out] first_reg should be passed in as the first region
         * index in the lattice, and will be incremented to give the first index
         * in the returned \ref PinMesh
         * \param[in] dir a \ref Direction used to disambiguate which \ref
         * PinMesh is desired when the \p p lies directly on a pin boundary. The
         * convention is to return the \ref PinMesh towards which \p dir points.
         *
         * These calls are chained from the CoreMesh -> Plane -> Lattice, with
         * each level in the geometrical hierarchy moving the point to the
         * appropriate local coordinates and offsetting the first_reg value.
         */
        const PinMesh* get_pinmesh( Point2 &p, int &first_reg,
                Direction dir=Direction() ) const;

        /**
         * \brief Return whether the current \ref Lattice and the passed \ref
         * Lattice reference are \ref Assembly compatible.
         *
         * \ref Assembly compatible means that the two \ref Lattice's can be
         * stacked atop each other and their pin boundaries will line up.
         */
        bool compatible( const Lattice &other ) const;

    private:
        size_t id_;
        unsigned nx_;
        unsigned ny_;
        size_t n_reg_;
        size_t n_xsreg_;
        real_t hx_;
        real_t hy_;
        VecF hx_vec_;
        VecF hy_vec_;
        VecF x_vec_;
        VecF y_vec_;

        // Array of pins in the lattice
        std::vector<Pin*> pins_;

        // Array of starting FSR indices for each pin in the lattice
        VecI first_reg_pin_;
    };

    typedef std::shared_ptr<Lattice> SP_Lattice_t;
    typedef std::shared_ptr<Lattice> UP_Lattice_t;

    std::map<int, UP_Lattice_t> ParseLattices( const pugi::xml_node &input,
            const std::map<int, UP_Pin_t> &pins );

}
