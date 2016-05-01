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

#include "pin_mesh_cyl.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "pugixml.hpp"

#include "constants.hpp"
#include "error.hpp"
#include "files.hpp"
#include "global_config.hpp"
#include "string_utils.hpp"

using std::string;
using std::stringstream;

using std::cout;
using std::endl;

namespace mocc {
    PinMesh_Cyl::PinMesh_Cyl(const pugi::xml_node &input):
        PinMesh( input ) {
        // Extract the radii and check for sanity
        {
            stringstream radiiIn(input.child("radii").child_value());
            while (!radiiIn.eof()) {
                mocc::real_t rad;
                radiiIn >> rad;
                xs_radii_.push_back(rad);

                if(radiiIn.fail()){
                    stringstream msg;
                    msg << "Ran into a problem reading radii for pin ID="
                        << id_;
                    throw EXCEPT(msg.str().c_str());
                }
            }
            if(!radiiIn.eof()){
                stringstream msg;
                msg << "Dangling data in radii for pin ID="
                    << id_;
            }
            // Make sure the radii are ordered
            for (unsigned int i=0; i<xs_radii_.size()-1; i++){
                if (xs_radii_[i] > xs_radii_[i+1]) {
                    // The radii are not ordered. Pitch a fit
                    stringstream msg;
                    msg << "Pin radii do not appear to be ordered for pin ID="
                        << id_;
                    throw EXCEPT(msg.str());
                }
            }

            // Make sure the last radius is smaller than a half-pitch
            if(xs_radii_.back() > pitch_x_*0.5){
                throw EXCEPT("Largest radius is too big!");
            }

            n_xsreg_ = xs_radii_.size() + 1;

        }

        // Read in the azimuthal subdivisions
        {
            sub_azi_ = VecI();
            stringstream inBuf(input.child("sub_azi").child_value());
            int azi;
            inBuf >> azi;
            if(inBuf.fail()){
                throw EXCEPT("Improper input to azimuthal subdivisions!");
            }
            sub_azi_.push_back(azi);

            // For now im only supporting one entry (same azi for all rings).
            if(sub_azi_.size() > 1){
                throw EXCEPT("Only supporting on azi type for now.");
            }

            // Make sure that the azimuthal division is even and <=8. One of
            // these days, ill solve the more general problem of any number of
            // azis.
            if( (sub_azi_[0]%2 != 0) | (sub_azi_[0] > 8) ) {
                throw EXCEPT("Only supporting even azimuthal subdivisions "
                        "<=8.");
            }
        }

        // Read in the radial subdivisions
        {
            sub_rad_ = explode_string<int>( 
                    input.child("sub_radii").child_value());

            // Make sure we have the same number of radial subdivs as rings
            if( (int)sub_rad_.size() != (n_xsreg_ - 1) ) {
                throw EXCEPT("Wrong number of radial subdivisions specified.");
            }
        }

        //
        // We should be done extracting information from the XML at this point
        //

        // Calculate actual mesh radii. They should have equal volume within each
        // XS ring.
        double rxsi = 0.0;
        double ri = 0.0;
        for( int ixs=0; ixs<n_xsreg_ - 1; ixs++ ) {
            double vn = (xs_radii_[ixs]*xs_radii_[ixs] - rxsi*rxsi) /
                        sub_rad_[ixs];

            for(int ir=0; ir<sub_rad_[ixs]; ir++){
                double r = sqrt(vn + ri*ri);
                radii_.push_back(r);
                ri = r;
            }
            rxsi = xs_radii_[ixs];
        }

        // Construct Circle objects corresponding to each mesh ring
        Point2 origin(0.0, 0.0);
        for( auto ri=radii_.begin(); ri!=radii_.end(); ++ri ) {
            circles_.push_back(Circle(origin, *ri));
        }

        // Construct Line objects corresponding to each azimuthal subdivision
        real_t h_pitch_x = 0.5*pitch_x_;
        real_t h_pitch_y = 0.5*pitch_y_;
        int n_azi = sub_azi_[0];
        Box pin_box( Point2(-h_pitch_x, -h_pitch_y),
                     Point2( h_pitch_x,  h_pitch_y) );
        real_t ang_sep = TWOPI/n_azi;
        for( int iazi=0; iazi<n_azi; iazi++ ) {
            Angle ang(iazi*ang_sep, HPI, 0.0);
            // Point at which the azimuthal subdivision intersects the bounding
            // box of the pin.
            Point2 p = pin_box.intersect(origin, ang);
            lines_.push_back( Line(origin, p) );
        }

        n_reg_ = (radii_.size()+1) * sub_azi_[0];


        // Determine FSR volumes
        real_t prev_r = 0.0;
        for( auto &ri: radii_) {
            real_t voli = PI*(ri*ri - prev_r*prev_r)/n_azi;
            for( int iazi=0; iazi<n_azi; iazi++ ) {
                vol_.push_back(voli);
            }
            prev_r = ri;
        }

        // Add the volumes from the outer annular region
        real_t large_r = radii_[radii_.size()-1];
        real_t v_outer = (pitch_x_*pitch_y_ - PI*large_r*large_r)/n_azi;
        for( int ia=0; ia<n_azi; ia++ ) {
            vol_.push_back(v_outer);
        }
        assert( (int)vol_.size() == n_reg_ );

        return;
    }

    PinMesh_Cyl::~PinMesh_Cyl() {
    }

    int PinMesh_Cyl::trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
            VecI &reg ) const {
        Line l(p1, p2);

        std::vector<Point2> ps;

        // Add the pin entry and exit points to the vector
        ps.push_back(p1);
        ps.push_back(p2);

        // Find intersections with the rings
        for( auto ci=circles_.begin(); ci!=circles_.end(); ci++ ) {
            Point2 p1;
            Point2 p2;
            int ret = Intersect( l, *ci, p1, p2 );
            if( ret == 2 ) {
                ps.push_back(p1);
                ps.push_back(p2);
            }
        }

        // Find intersections with the azimuthal subdivisions
        for( auto &li: lines_ ) {
            Point2 p;
            int ret = Intersect( li, l, p );
            if( ret == 1 ) {
                ps.push_back(p);
            }
        }

        // Sort the intersection points and remove duplicates
        std::sort(ps.begin(), ps.end());
        ps.erase( std::unique(ps.begin(), ps.end()), ps.end() );

        // Determine segment lengths and region indices
        for( unsigned int ip=1; ip<ps.size(); ip++ ) {
            s.push_back( ps[ip].distance(ps[ip-1]) );
            unsigned int local_reg =
                this->find_reg( Midpoint(ps[ip], ps[ip-1]) );
            reg.push_back( local_reg + first_reg );
        }
        return ps.size()-1;
    }

    /**
     * For now, indexing in the cylindrical pins goes from the inside radius
     * out, and from the positive x axis around azimuthally counter-clockwise.
     * At some point, I might look into other indexing schemes to try and
     * achieve better locality and cache performance, but for now KISS.
     */
    int PinMesh_Cyl::find_reg( Point2 p ) const {
        // Test that the point is inside the pin mesh
        if( (fabs(p.x) > 0.5*pitch_x_) | (fabs(p.y) > 0.5*pitch_y_) ) {
            return -1;
        }

        // Find the radial division of the point
        real_t r = sqrt( p.x*p.x+p.y*p.y );
        unsigned int ir = 0;
        for( ir=0; ir<radii_.size(); ir++ ) {
            if( r < radii_[ir]) {
                break;
            }
        }

        // This is only a little tricky; if the loop above runs through it means
        // that the point is outside the largest ring, and therefore in the
        // annular region outside the pin. Conveniently, ir will be the proper
        // index corresponding to that region, so we can go ahead and use it.

        // Find the azimuthal subdivision that the point is in.
        real_t azi = p.alpha();
        int ia = azi/(TWOPI/sub_azi_[0]);
        int ireg = ir*sub_azi_[0] + ia;

        assert( (0 <= ireg) & (ireg < n_reg_ ) );

        return ireg;
    }

    void PinMesh_Cyl::print( std::ostream &os ) const {
        PinMesh::print( os );
        os << std::endl;
        os << "Type: Cylindrical";
        return;
    }

    std::string PinMesh_Cyl::draw() const {
        std::stringstream buf;

        buf << "ctx.move_to(0, 0)" << std::endl;
        for( auto c: circles_ ) {
            buf << "ctx.arc(" << c.c.x << ", "
                              << c.c.y << ", "
                              << c.r << ", 0, twopi)" << std::endl;;
        }

        for( auto l: lines_ ) {
            buf << "ctx.move_to(" << l.p1.x << ", "
                                  << l.p1.y << ")" << std::endl;
            buf << "ctx.line_to(" << l.p2.x << ", "
                                  << l.p2.y << ")" << std::endl;
            buf << "ctx.close_path()" << std::endl;
        }
        buf << "ctx.stroke()";
        return buf.str();
    }
}
