#include "pin_mesh_rect.hpp"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>

#include "error.hpp"

namespace mocc {
    PinMesh_Rect::PinMesh_Rect( const pugi::xml_node &input ):
        PinMesh( input )
    {
        // Parse the number of x and y divisions
        int ndiv_x = input.child("sub_x").text().as_int(0);
        if (ndiv_x < 1) {
            throw EXCEPT("Failed to read valid number of X divisions in rect "
                    "pin mesh.");
        }
        int ndiv_y = input.child("sub_y").text().as_int(0);
        if (ndiv_y < 1) {
            throw EXCEPT("Failed to read valid number of Y divisions in rect "
                    "pin mesh.");
        }

        n_xsreg_ = ndiv_x * ndiv_y;
        n_reg_   = ndiv_x * ndiv_y;

        real_t dx = pitch_x_/ndiv_x;
        real_t dy = pitch_y_/ndiv_y;

        real_t h_pitch_x = 0.5*pitch_x_;
        real_t h_pitch_y = 0.5*pitch_y_;

        for (int i=1; i<ndiv_x; i++) {
            hx_.push_back( i*dx-h_pitch_x );
        }
        for (int i=1; i<ndiv_y; i++) {
            hy_.push_back( i*dy-h_pitch_y );
        }

        // Form lines representing the mesh boundaries

        for ( auto &xi: hx_ ) {
            lines_.push_back( Line( Point2(xi, -h_pitch_y),
                                    Point2(xi,  h_pitch_y) ) );
        }
        for ( auto &yi: hy_ ) {
            lines_.push_back( Line( Point2(-h_pitch_x, yi),
                                    Point2( h_pitch_x, yi) ) );
        }

        // Determine FSR volumes
        vol_ = VecF(n_reg_, dx*dy);

        return;
    }

    int PinMesh_Rect::trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
            VecI &reg ) const {

        // Make a vector to store the collision points and add the start and end
        // points to it
        std::vector<Point2> ps;
        ps.push_back(p1);
        ps.push_back(p2);

        // Create a line object for the input points to test for collisions with
        Line l(p1, p2);

        // Test for collisions with all of the lines in the mesh
        for (auto &li: lines_) {
            Point2 p;
            int ret = Intersect( li, l, p );
            if( ret == 1 ) {
                ps.push_back(p);
            }
        }

        // Sort the points
        std::sort( ps.begin(), ps.end() );
        ps.erase( std::unique(ps.begin(), ps.end()), ps.end() );

        // Determine segment lengths and region indices
        for( unsigned int ip=1; ip<ps.size(); ip++ ) {
            s.push_back( ps[ip].distance(ps[ip-1]) );
            reg.push_back( this->find_reg( Midpoint(ps[ip], ps[ip-1]) )
                    + first_reg );
        }

        return ps.size()-1;
    }

    // Return an integer containing the pin-local FSR index in which the passed
    // Point resides.
    //
    // In the rectangular mesh, the indices are ordered naturally. The first
    // region is in the lower left, the last in the upper right, proceeding
    // first in the x direction, then in the y. nothing too fancy.
    int PinMesh_Rect::find_reg( Point2 p ) const {
        // Make sure the point is inside the pin
        if( fabs(p.x) > 0.5*pitch_x_ ) {
            return -1;
        }
        if( fabs(p.y) > 0.5*pitch_y_ ) {
            return -1;
        }

        int ix;
        for( ix=0; ix<(int)hx_.size(); ix++ ) {
            if( hx_[ix] > p.x ) {
                break;
            }
        }
        int iy;
        for( iy=0; iy<(int)hy_.size(); iy++ ) {
            if( hy_[iy] > p.y ) {
                break;
            }
        }

        int ireg = (hx_.size()+1)*iy + ix;
        assert(ireg < n_reg_);
        return ireg;
    }

    void PinMesh_Rect::print( std::ostream &os ) const {
        PinMesh::print( os );
        os << std::endl;
        os << "Type: Rectangular";
        return;
    }

    std::string PinMesh_Rect::draw() const {
        std::stringstream buf;

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
