#include "ray.hpp"


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

        mesh.trace( ps, cm_start_fw_, cm_start_bw_, s_fw, s_bw );

        // Now we have a list of points of intersections with all of the pin
        // boundaries. 
        // Start by finding the first coarse surface crossing
        Point2 pin_p = Midpoint(ps[0], ps[1]);
        int cell = mesh.coarse_cell_point(pin_p);
        int surf[2];
        int nsurf = mesh.coarse_surf_point( ps[0], cell, surf );
        for( int i=0; i<nsurf; i++ ) {
            cm_surf_.push_back(surf[i]);
        }

        // Loop over them and trace individual pin meshes.
        Point2 p_prev = ps[0];
        for( auto pi=ps.begin()+1; pi!=ps.end(); ++pi ) {
            // Use the midpoint of the pin entry and exit points to locate the
            // pin.
            pin_p = Midpoint(*pi, p_prev);

            int first_reg = 0;
            const PinMeshTuple pmt = mesh.get_pinmesh(pin_p, iz, first_reg);

            int nseg = pmt.pm->trace(p_prev-pin_p, *pi-pin_p, first_reg,
                    seg_len_, seg_index_);


            // Figure out coarse mesh info.
            cm_nseg_.push_back(nseg);
            cm_cell_.push_back( mesh.coarse_cell(pmt.position) );
            int surf[2];
            int nsurf = mesh.coarse_surf_point( *pi, cm_cell_.back(), 
                    surf );
            for( int i=0; i<nsurf; i++ ) {
//cout << surf[i] << endl;
                cm_surf_.push_back(surf[i]);
            }
            
            p_prev = *pi;
        }

        nseg_ = seg_len_.size();

        bc_[0] = bc1;
        bc_[1] = bc2;

        return;
    }

}
