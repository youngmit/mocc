#include "ray.hpp"

#include <cassert>

using std::cout;
using std::endl;

// Assuming that p1 is the "origin" return the quadrant of the angle formed by
// p1. Since we assume that p1 is below p2 in y, only octants 1 or 2 can be
// returned.
inline int get_octant( mocc::Point2 p1, mocc::Point2 p2 ) {
    assert(p2.y > p1.y);
    return ( p2.x > p1.x ) ? 1 : 2;
}


namespace mocc {
    /** 
     * \param p1 the starting point of the Ray.
     * \param p2 the ending point of the Ray.
     * \param bc1 the boundary condition index corresponding to p1.
     * \param bc2 the boundary condition index corresponding to p2.
     * \param iz the index of the geometry to trace. This can be any Z index
     * that contains the geometrically-unique place that we are generating ray
     * data for.
     * \param mesh a reference to the CoreMesh to trace.
     *
     * A Ray is defined by two Point2 structs, specifying the beginning and end
     * of a ray, on the boundary of the problem. Given these two points, all of
     * the segments of the ray are determined by first finding intersections
     * with the pin cell edges (using CoreMesh::trace()), then the internal
     * surface crossings for each pin (using PinMesh::trace()).
     */
    Ray::Ray( Point2 p1, Point2 p2, size_t bc1, size_t bc2, int iz, 
            const CoreMesh &mesh ) {
        std::vector<Point2> ps;

        ps.push_back(p1);
        ps.push_back(p2);

        std::vector<Surface> s_fw;
        std::vector<Surface> s_bw;
        VecI cm_nseg;

        mesh.trace( ps );

        // Determine the octant of the angle of the ray
        int octant_fw = get_octant( p1, p2 );
        int octant_bw = (octant_fw == 1) ? 3 : 4;
        cm_start_fw_ = mesh.coarse_boundary_cell( p1, octant_fw );
        cm_start_bw_ = mesh.coarse_boundary_cell( p2, octant_bw );
        
        size_t ns = 0;
        
        // Loop over them and trace individual pin meshes.
        Point2 p_prev = ps[0];

        // Always start with zero segments at the beginning, since the first
        // surface normal is "upwind"
        cm_nseg.push_back( 0 );

        // Get the starting CM ray info for FW, ending for BW
        Surface s[2];
        ns = mesh.coarse_norm_point( p_prev, octant_fw, s );
        for( int i=0; i<ns; i++ ) {
            s_fw.push_back(s[i]);
        }
        if( ns == 2 ) {
            cm_nseg.push_back( 0 );
        }
        ns = mesh.coarse_norm_point( p_prev, octant_bw, s );
        for( int i=0; i<ns; i++ ) {
            s_bw.push_back(s[i]);
        }


        for( auto pi=ps.begin()+1; pi!=ps.end(); ++pi ) {
            // Use the midpoint of the pin entry and exit points to locate the
            // pin.
            auto pin_p = Midpoint(*pi, p_prev);

            int first_reg = 0;
            const PinMeshTuple pmt = mesh.get_pinmesh(pin_p, iz, first_reg);

            int nseg = pmt.pm->trace(p_prev-pin_p, *pi-pin_p, first_reg,
                    seg_len_, seg_index_);


            ns = mesh.coarse_norm_point( *pi, octant_fw, s );
            for( int i=0; i<ns; i++ ) {
                s_fw.push_back( s[i] );
            }
            cm_nseg.push_back( nseg );
            if( ns == 2 ) {
                cm_nseg.push_back( 0 );
            }
            ns = mesh.coarse_norm_point( *pi, octant_bw, s );
            for( int i=0; i<ns; i++ ) {
                s_bw.push_back( s[i] );
            }
            
            p_prev = *pi;
        }

        // Flip the order of the BW surfaces
        std::reverse( s_bw.begin(), s_bw.end() );

        /**
         * \todo there is room to do a couple clevers here, maybe only store one
         * value for nseg. Should come back to it when i get a chance.
         */

        // Here ns is the offset by which we should access the number of
        // segments in the backward direction from the list of segments in the
        // forward direction. Normally, it should be size-1, or the last
        // element. If the ray crosses a corner at the end in the FW direction,
        // there is zero that we need to jump over.
        ns = cm_nseg.size();
        if( cm_nseg[cm_nseg.size()-1] == 0 ) {
            ns -= 1;
        }
        for( size_t i=0; i<cm_nseg.size(); i++ ) {
            RayCoarseData rcd;
            rcd.fw = s_fw[i];
            rcd.bw = s_bw[i];
            rcd.nseg_fw = cm_nseg[i];
            rcd.nseg_bw = cm_nseg[ (ns-i)%cm_nseg.size() ];
            cm_data_.push_back(rcd);
        }
        /*
        {
            RayCoarseData rcd;
            rcd.fw = s_fw.back();
            rcd.bw = s_bw.back();
            rcd.nseg_fw = 0;
            rcd.nseg_bw = 0;
            cm_data_.push_back(rcd);
        }
        */

        nseg_ = seg_len_.size();

        bc_[0] = bc1;
        bc_[1] = bc2;
        return;
    }

}
