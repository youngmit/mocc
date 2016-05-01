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

#include <string>

#include "global_config.hpp"
#include "geom.hpp"
#include "position.hpp"
#include "pugifwd.hpp"

namespace mocc {
    /**
     * \ref PinMesh is a virtual class, which provides methods for performing
     * ray tracing and accessing data in common between all types of pin mesh,
     * such as region volumes, x and y pitch, etc.
    */
    class PinMesh {
    public:
        PinMesh( const pugi::xml_node &input );

        virtual ~PinMesh() {
        }

        int id() const {
            return id_;
        }

        int n_reg() const {
            return n_reg_;
        }

        int n_xsreg() const {
            return n_xsreg_;
        }

        real_t pitch_x() const {
            return pitch_x_;
        }

        real_t pitch_y() const {
            return pitch_y_;
        }

        real_t vol() const {
            return pitch_x_ * pitch_y_;
        }

        const VecF& vols() const {
            return vol_;
        }

        /**
         * \param[in] p1 the entry point of the ray.
         * \param[in] p2 the exit point of the ray.
         * \param[in] first_reg the index of the first FSR in the pin mesh.
         * \param[in,out] s the vector of ray segment lengths to be modified.
         * \param[in,out] reg the vector of ray segment FSR indices to be
         * modified.
         *
         * \returns The number of segments that pass through pin geometry
         * (useful for CMFD data).
         *
         * Given an entry and exit point, which should be on the boundary of the
         * pin (in pin-local coordinates), and the index of the first FSR in the
         * pin, append values to the vectors of segment length and corresponding
         * region index.
         *
         *
         * The segment lengths are uncorrected, which is to say that they are
         * the true lengths of the rays as they pass through the mesh.
         * Therefore, summing the volume of the segments in each FSR is not
         * guaranteed to return the correct FSR volume. Make sure to correct for
         * this after tracing all of the rays in a given angle.
        */
        virtual int trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
                VecI &reg ) const =0;

        /**
         * Find the pin-local region index corresponding to the point provided.
         *
         * \param[in] p the \ref Point2 for which to find the FSR index
         *
         * Given a point in pin-local coordinates, return the mesh region index
         * in which the point resides
        */
        virtual int find_reg( Point2 p ) const =0;

        /**
         * Return the number of flat source regions corresponding to an XS
         * region (indexed pin-locally).
        */
        virtual size_t n_fsrs( unsigned int xsreg ) const =0;

        /**
         * \brief Insert hte \ref PinMesh into the passed ostream.
         *
         * This is useful to be able to use \c operator<< on a polymorphic \ref
         * PinMesh and have it work as expected.
         */
        virtual void print( std::ostream &os ) const;

        /**
         * \brief Return a string containing PyCairo commands to draw the \ref
         * PinMesh.
         */
        virtual std::string draw() const = 0;

        /**
         * \brief Provide stream insertion support.
         *
         * Calls the \ref PinMesh::print() routine, which is virtual, allowing
         * for polymorphism.
         */
        friend std::ostream& operator<<( std::ostream &os, const PinMesh& pm );

    protected:
        int id_;
        int n_reg_;
        int n_xsreg_;
        real_t pitch_x_;
        real_t pitch_y_;
        VecF vol_;
    };

    /**
     * This is a simple struct that contains a const pointer to a \ref PinMesh,
     * along with a Position describing its location. This is essentially a
     * useful tuple for returning both values from a lookup function (see
     * \ref CoreMesh::get_pinmesh() and \ref Plane::get_pinmesh()).
    */
    struct PinMeshTuple {
        PinMeshTuple( Position pos, const PinMesh *pm ):
            position( pos ),
            pm( pm )
        { }
        Position position;
        const PinMesh *pm;
    };
}
