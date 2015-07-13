#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "global_config.hpp"
#include "geom.hpp"

namespace mocc {
    class PinMesh {
    public:
        PinMesh( const pugi::xml_node &input );

        virtual ~PinMesh() {
        }
    	
        unsigned int id() const {
    		return id_;
    	}

        unsigned int n_reg() const {
            return n_reg_;
        }

        unsigned int n_xsreg() const {
            return n_xsreg_;
        }

        float_t pitch_x() const {
            return pitch_x_;
        }

        float_t pitch_y() const {
            return pitch_y_;
        }

        const VecF& vol() const {
            return vol_;
        }

        // Given an entry and exit point [, whichshould be] on the boundary of
        // the pin, in pin-local coordinates, and the first region index, append
        // values to the vectors of segment length and region index. 
        //
        // Returns the number of segments that pass through pin geometry (for
        // CMFD data).
        //
        // The segment lengths are uncorrected, which is to say that they are
        // the true lengths of the rays as they pass through the mesh.
        // Therefore, summing the volume of the segments in each FSR is not
        // guaranteed to return the correct FSR volume. Make sure to correct for
        // this after stracing all of the rays in a given angle.
        virtual int trace( Point2 p1, Point2 p2, int first_reg, VecF &s, 
                VecI &reg ) const =0;

        // Given a point in pin-local coordinates, return the mesh region index
        // in which the point resides
        virtual int find_reg( Point2 p ) const =0;

    protected:
    	unsigned int id_;
        unsigned int n_reg_;
        unsigned int n_xsreg_;
    	float_t pitch_x_;
    	float_t pitch_y_;
        VecF vol_;
    };
}
