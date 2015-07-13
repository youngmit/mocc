#include "ray_data.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

#include "error.hpp"

using std::cout;
using std::endl;

namespace mocc {

    // Construct a ray from a starting point and angle
    Ray::Ray( Point2 p1, Point2 p2, int iz, const CoreMesh &mesh ) {
        std::vector<Point2> ps;

        ps.push_back(p1);
        ps.push_back(p2);

        mesh.trace( ps );

        // Now we have a list of points of intersections with all of the pins.
        // Loop over them and trace individual pin meshes.
        Point2 p_prev = ps[0];
        for (auto pi=ps.begin()+1; pi!=ps.end(); ++pi) {
            Point2 pin_p = Midpoint(*pi, p_prev);

            int first_reg = 0;
            const PinMesh* pm = mesh.get_pinmesh(pin_p, iz, first_reg);

            pm->trace(p_prev-pin_p, *pi-pin_p, first_reg, seg_len_, seg_index_);
            p_prev = *pi;
        }

        nseg_ = seg_len_.size();

        return;
    }


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

            float_t new_alpha = atan( hy*Nx/(hx*Ny) );
            ang = ModifyAlpha( ang, new_alpha );

            ang_quad_.modify_angle( iang, ang );
            float_t space = cos(ang_it->alpha) * hy/Ny;
cout << "spacing: " << space << endl;
            spacing_.push_back( space );

            iang++;
        }

cout << "Modularized Angular quadrature " << endl;
cout << ang_quad_ << endl;


        // Trace rays
        Box core_box = Box(Point2(0.0, 0.0), Point2(hx, hy));
        // loop over the planes of unique geometry
        for ( unsigned int iplane=0; iplane<n_planes_; iplane++ ) {
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
                int Nx = Nx_[iang%ang_quad_.ndir_oct()];
                int Ny = Ny_[iang%ang_quad_.ndir_oct()];
                float_t space = spacing_[iang%ang_quad_.ndir_oct()];
                float_t space_x = fabs( space/sin(ang->alpha) );
                float_t space_y = fabs( space/cos(ang->alpha) );

                std::vector<Ray> rays;
                // Handle rays entering on the x-normal face
                for ( int iray=0; iray<Ny; iray++ ) {
                    Point2 p1;
                    if( ang->ox > 0.0 ) {
                        // We are in octant 1, enter from the left/west
                        p1.x = 0.0;
                    } else {
                        // We are in octant 2, enter from the right/east
                        p1.x = hx;
                    }
                    p1.y = (0.5 + iray)*space_y;
                    Point2 p2 = core_box.intersect(p1, *ang);
                    rays.emplace_back(Ray(p1, p2, iplane, mesh));
                }
                for ( int iray=0; iray<Nx; iray++ ) { 
                    Point2 p1;
                    p1.x = (0.5 + iray)*space_x;
                    p1.y = 0.0;
                    Point2 p2 = core_box.intersect(p1, *ang);
                    rays.emplace_back(Ray(p1, p2, iplane, mesh));
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
                }

                // Move the stack of rays into the vector of angular ray sets.
                angle_rays.push_back(std::move(rays));
                ++iang;
            }
            // Move the angular ray set to the vector of planar ray sets.
            rays_.push_back(std::move(angle_rays));
        }

        // Adjust ray lengths to correct FSR volume. Use an angle integral to do
        // so.
        for( unsigned int iplane=0; iplane<n_planes_; iplane++ ) {
            const VecF& true_vol = mesh.plane(iplane).vol();
            VecF fsr_vol(mesh.plane(iplane).n_reg(), 0.0);
            int iang=0;
            for( auto ang = ang_quad_.octant(1); 
                 ang!=ang_quad_.octant(3);
                 ++ang ) {
                std::vector<Ray>& rays = rays_[iplane][iang];
                float_t space = spacing_[iang%ang_quad_.ndir_oct()];
                float_t wgt = ang->weight*0.5;

                for( auto &ray: rays ) {
                    for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        unsigned int ireg = ray.seg_index(iseg);
                        fsr_vol[ireg] += ray.seg_len(iseg) * space * wgt ;
                    }
                }
                ++iang;
            }
            // Convert fsr_vol into a correction factor
            for( unsigned int ireg=0; ireg<mesh.plane(iplane).n_reg(); ireg++ ) {
                fsr_vol[ireg] = true_vol[ireg]/fsr_vol[ireg];
            }
            // Correct ray lengths to enforce proper FSR volumes
            iang = 0;
            for( auto ang = ang_quad_.octant(1); 
                 ang!=ang_quad_.octant(3);
                 ++ang ) {
                std::vector<Ray>& rays = rays_[iplane][iang];
                for( auto &ray: rays ){
                    for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        unsigned int ireg = ray.seg_index(iseg);
                        ray.seg_len(iseg) = ray.seg_len(iseg)*fsr_vol[ireg];
                    }
                }
                ++iang;
            }
        } // Volume correction

    }


}
