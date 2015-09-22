#include "ray_data.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

#include "error.hpp"
#include "constants.hpp"

using std::cout;
using std::endl;
using std::cin;

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

        mesh.trace( ps );

        // Now we have a list of points of intersections with all of the pin
        // boundaries. 
        // Start by finding the first coarse surface crossing
        Point2 pin_p = Midpoint(ps[0], ps[1]);
        int cell = mesh.coarse_cell_point(pin_p);
        int surf[2];
        int nsurf = mesh.coarse_surf_point( ps[0], cell, surf );
        if( nsurf == 1 ) {
            cm_surf_.push_back(surf[0]);
        } else {
            throw EXCEPT("Boundary point should be unique.");
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
                cm_surf_.push_back(surf[i]);
            }
            
            p_prev = *pi;
        }

        nseg_ = seg_len_.size();

        bc_[0] = bc1;
        bc_[1] = bc2;

        return;
    }

    /**
     * \brief Construct a RayData object using a \<rays\> XML tag, a desired
     * AngularQuadrature, and CoreMesh.
     *
     * For now, the angular quadrature is duplicated before performing
     * modularization, which only mutates the RayData quadrature.
     *
     * \todo At some point, it is necessary to get access to the modularized
     * quadrature. There are a couple options for this:
     * - Define the modularization procedure on the AngularQuadrature itself,
     *   and perform modularization on the quadrature before we even get to this
     *   point. I don't like this, becuase it doesnt make sense conceptually to
     *   live on the quadrature itself.
     * - Provide a means to return the modularized quadrature after construction
     *   of the RayData. A little better, but perhaps cumbersome.
     * - Allow the passed AngularQuadrature reference to be non-const and mutate
     *   it directly.
     *
     * Construction performs the following steps:
     * -# Parse input from the XML
     * -# Modularize the angular quadrature and determine ray spacing for each
     *  angle
     * -# Construct Ray objects for each geometrically-unique plane and angle
     * -# Correct the ray segment lengths to preserve FSR volumes
     *
    */
    RayData::RayData( const pugi::xml_node &input, 
            const AngularQuadrature &ang_quad,
            const CoreMesh &mesh ):
    ang_quad_(ang_quad)
    {
        // Make sure we have reasonable input
        if ( input.empty() ) {
            Error("No input privided for ray spacing.");
        }

        // Get the optimal ray spacing
        float_t opt_spacing = input.attribute("spacing").as_float(-1.0);
        if( opt_spacing <= 0.0 ) {
            Error("Failed to read valid ray spacing.");
        }

        // Store some necessary stuff from the CoreMesh
        n_planes_ = mesh.n_unique_planes();

        // Figure out modular angles and spacings
        float_t hx = mesh.hx();
        float_t hy = mesh.hy();

cout << "Original Angular quadrature " << endl;
cout << ang_quad_ << endl;

        int iang = 0;
        for (auto ang_it = ang_quad_.octant(1); 
                ang_it != ang_quad_.octant(2); ++ang_it) {

            Angle ang = *ang_it;

            int Nx = ceil( hx/opt_spacing*fabs( sin( ang.alpha ) ) );
            int Ny = ceil( hy/opt_spacing*fabs( cos( ang.alpha ) ) );
            Nx += Nx%2+1;
            Ny += Ny%2+1; 

            Nx_.push_back(Nx);
            Ny_.push_back(Ny);
            Nrays_.push_back(Nx+Ny);

            float_t new_alpha = atan( hy*Nx/(hx*Ny) );
            ang = ModifyAlpha( ang, new_alpha );

            ang_quad_.modify_angle( iang, ang );
            float_t space = cos(ang_it->alpha) * hy/Ny;
            spacing_.push_back( space );

            iang++;
        }
        
        // push more Nx, Ny, N, space onto their respective vectors, so we dont
        // have to wory about %'ing by ndir_oct
        for( iang=0; iang<ang_quad_.ndir_oct()*3; iang++ ) {
            Nx_.push_back(Nx_[iang]);
            Ny_.push_back(Ny_[iang]);
            Nrays_.push_back(Nrays_[iang]);
            spacing_.push_back(spacing_[iang]);
        }

cout << "Modularized Angular quadrature " << endl;
cout << ang_quad_ << endl;

        // Trace rays
        Box core_box = Box(Point2(0.0, 0.0), Point2(hx, hy));
        max_seg_ = 0;
        // loop over the planes of unique geometry
        for ( size_t iplane=0; iplane<n_planes_; iplane++ ) {
            // generate rays for each angle in octants 1 and 2
            int iang = 0;
            std::vector< std::vector<Ray> > angle_rays;
            int nreg_plane = mesh.plane(iplane).n_reg();
            for ( auto ang = ang_quad_.octant(1);
                    ang != ang_quad_.octant(3); ++ang ) {
                VecI nrayfsr( nreg_plane, 0 );
                // We only define Nx, Ny and spacing above for the first octant,
                // so mod the angle index with the number of angles per octant
                // to cast it into the first octant.
                int Nx = Nx_[iang];
                int Ny = Ny_[iang];
                int bc1 = 0;
                int bc2 = 0;
                float_t space = spacing_[iang];
                float_t space_x = fabs( space/sin(ang->alpha) );
                float_t space_y = fabs( space/cos(ang->alpha) );

                std::vector<Ray> rays;
                // Handle rays entering on the x-normal faces
                for ( int iray=0; iray<Ny; iray++ ) {
                    Point2 p1;
                    bc1 = iray;
                    if( ang->ox > 0.0 ) {
                        // We are in octant 1, enter from the left/west
                        p1.x = 0.0;
                    } else {
                        // We are in octant 2, enter from the right/east
                        p1.x = hx;
                    }
                    p1.y = (0.5 + iray)*space_y;
                    Point2 p2 = core_box.intersect(p1, *ang);
                    if( fp_equiv_abs(p2.x, hx) ) {
                        // BC is on the right/east boundary of the core
                        bc2 = p2.y/space_y;
                    } else if ( fp_equiv_abs(p2.y, hy) ) {
                        // BC is on the top/north boundary of the core
                        bc2 = p2.x/space_x + Ny;
                    } else if ( fp_equiv_abs(p2.x, 0.0) ) {
                        // BC is on the left/west boundary of the core
                        bc2 = p2.y/space_y;
                    } else {
                        Error("Something has gone horribly wrong in the ray "
                              "trace.");
                    }
                    rays.emplace_back(Ray(p1, p2, bc1, bc2, iplane, mesh));
                    max_seg_ = std::max( rays.back().nseg(), max_seg_ );
                }

                // Handle rays entering on the y-normal face
                for ( int iray=0; iray<Nx; iray++ ) { 
                    Point2 p1;
                    p1.x = (0.5 + iray)*space_x;
                    p1.y = 0.0;
                    Point2 p2 = core_box.intersect(p1, *ang);
                    bc1 = iray+Ny;
                    if( fp_equiv_abs(p2.x, hx) ) {
                        // BC is on the right/east boundary of the core
                        bc2 = p2.y/space_y;
                    } else if( fp_equiv_abs(p2.y, hy) ) {
                        // BC is on the top/north boundary of the core
                        bc2 = p2.x/space_x + Ny;
                    } else if( fp_equiv_abs(p2.x, 0.0) ) {
                        // BC is on the left/west boundary of the core
                        bc2 = p2.y/space_y;
                    } else {
                        Error("Something has gone horribly wrong in the ray "
                              "trace.");
                    }
                    rays.emplace_back(Ray(p1, p2, bc1, bc2, iplane, mesh));
                    max_seg_ = std::max( rays.back().nseg(), max_seg_ );
                }

                // Count number of ray crossings in each FSR
                for( auto &r: rays ) {
                    for ( auto &i: r.seg_index() ) {
                        nrayfsr[i]++;
                    }
                }

                // Make sure that there is at least one ray in every FSR. Give a
                // warning if not.
                if ( std::any_of(nrayfsr.begin(), nrayfsr.end(), 
                        [](int i){return i==0;}) ) {
                    Warn("No rays passed through at least one FSR. Try finer "
                            "ray spacing or larger regions.");
                    for( size_t ifsr=0; ifsr<nrayfsr.size(); ifsr++) {
                        cout << ifsr << " " << nrayfsr[ifsr] << endl;
                    }
                }

                // Move the stack of rays into the vector of angular ray sets.
                angle_rays.push_back(std::move(rays));
                ++iang;
            } // Angle loop
            // Move the angular ray set to the vector of planar ray sets.
            rays_.push_back(std::move(angle_rays));
        }

        // Adjust ray lengths to correct FSR volume. Use an angle integral to do
        // so.
        this->correct_volume( mesh, FLAT );
    }

    void RayData::correct_volume( const CoreMesh& mesh, VolumeCorrection type ) 
    {
        switch(type) {
            case FLAT:
                for( size_t iplane=0; iplane<n_planes_; iplane++ ) {
                    const VecF& true_vol = mesh.plane(iplane).vols();
                    int iang=0;
                    for ( auto ang = ang_quad_.octant(1); 
                            ang!=ang_quad_.octant(3); ++ang ) 
                    {
                        VecF fsr_vol(mesh.plane(iplane).n_reg(), 0.0);
                        std::vector<Ray>& rays = rays_[iplane][iang];
                        float_t space = spacing_[iang];
                        
                        for( auto &ray: rays ) {
                            for( size_t iseg=0; iseg<ray.nseg(); iseg++ )
                            {
                                size_t ireg = ray.seg_index(iseg);
                                fsr_vol[ireg] += ray.seg_len(iseg) * space;
                            }
                        }

                        // Correct
                        for( auto &ray: rays ) {
                            for( size_t iseg=0; iseg<ray.nseg(); iseg++ )
                            {
                                size_t ireg = ray.seg_index(iseg);
                                ray.seg_len(iseg) = ray.seg_len(iseg) *
                                    true_vol[ireg]/fsr_vol[ireg];
                            }
                        }
                        iang++;
                    }
                }

                break;
            case ANGLE:
                for( size_t iplane=0; iplane<n_planes_; iplane++ ) {
                    const VecF& true_vol = mesh.plane(iplane).vols();
                    VecF fsr_vol(mesh.plane(iplane).n_reg(), 0.0);
                    int iang=0;
                    for( auto ang = ang_quad_.octant(1); 
                         ang!=ang_quad_.octant(3); ++ang )
                    {
                        std::vector<Ray>& rays = rays_[iplane][iang];
                        float_t space = spacing_[iang];
                        float_t wgt = ang->weight*0.5;

                        for( auto &ray: rays ) {
                            for( size_t iseg=0; iseg<ray.nseg(); iseg++ )
                            { 
                                size_t ireg = ray.seg_index(iseg);
                                fsr_vol[ireg] += ray.seg_len(iseg) * space * 
                                    wgt;
                            }
                        }
                        ++iang;
                    }
                    // Convert fsr_vol into a correction factor
                    for( size_t ireg=0; ireg<mesh.plane(iplane).n_reg();
                            ireg++ ) 
                    {
                        fsr_vol[ireg] = true_vol[ireg]/fsr_vol[ireg];
                    }
                    // Correct ray lengths to enforce proper FSR volumes
                    iang = 0;
                    for( auto ang = ang_quad_.octant(1); 
                         ang!=ang_quad_.octant(3);
                         ++ang ) {
                        std::vector<Ray>& rays = rays_[iplane][iang];
                        for( auto &ray: rays ){
                            for( size_t iseg=0; iseg<ray.nseg(); iseg++ )
                            { 
                                size_t ireg = ray.seg_index(iseg);
                                ray.seg_len(iseg) = ray.seg_len(iseg) * 
                                    fsr_vol[ireg];
                            }
                        }
                        ++iang;
                    }
                } // Volume correction
                break;
        
        }
    }
}
