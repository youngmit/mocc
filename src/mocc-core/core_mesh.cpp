#include "core_mesh.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

#include "error.hpp"
#include "files.hpp"
#include "string_utils.hpp"


using std::cout;
using std::endl;
using std::stringstream;


namespace mocc {
    /**
     * \todo Simplify the CoreMesh construction by offloading stuff to the other
     * structure (Lattice, Core, etc.).
     */
    CoreMesh::CoreMesh(pugi::xml_node &input) {
        // Parse meshes
        pin_meshes_ = ParsePinMeshes( input );

        // Parse Material Library
        mat_lib_ = MaterialLib( input.child("material_lib") );

        // Parse pins
        pins_ = ParsePins( input, pin_meshes_ ); 

        // Parse lattices
        lattices_ = ParseLattices( input, pins_ );

        // Parse assemblies
        assemblies_ = ParseAssemblies( input, lattices_ );

        // Parse core
        core_ = ParseCore( input, assemblies_ );

        nx_ = core_.npin_x();
        ny_ = core_.npin_y();
        nz_ = core_.nz();
        n_surf_plane_ = (nx_+1)*ny_ + (ny_+1)*nx_ + nx_*ny_;
        nasy_ = core_.nasy();
        bc_ = core_.boundary();

        // Calculate the total core dimensions
        hx_ = 0.0;
        for ( int ix=0; ix<core_.nx(); ix++ ) {
            hx_ += core_.at(ix, 0).hx();
        }
        hy_ = 0.0;
        for ( int iy=0; iy<core_.ny(); iy++ ) {
            hy_ += core_.at(0, iy).hy();
        }

        // Determine the set of geometricaly-unique axial planes
        std::vector< VecI > unique;
        VecI plane_pins;
        int plane_reg = 0;
        for ( unsigned int iz=0; iz<nz_; iz++) {
            first_reg_plane_.push_back(plane_reg);
            // Form a list of all pin meshes in the core plane iz
            for ( unsigned int iasy=0; iasy<nasy_; iasy++ ) {
                const Assembly& asy = core_.at(iasy);
                for ( auto &pin: asy[iz] ) {
                    plane_pins.push_back(pin->mesh_id());
                    core_pins_.push_back(pin);
                    plane_reg += pin->n_reg();
                }
            }
            // Check against current list of unique planes
            int match_plane = -1;
            for( unsigned int iiz=0; iiz<unique.size(); iiz++ ) {
                for( unsigned int ip=0; ip<plane_pins.size(); ip++ ) {
                    if( plane_pins[ip] != unique[iiz][ip] ) {
                        // we dont have a match
                        break;
                    }
                    if ( ip == plane_pins.size()-1 ) {
                        // Looks like all of the pins matched!
                        match_plane = iiz;
                    }
                }
                if ( match_plane != -1 ) {
                    break;
                }
            }
            if ( match_plane == -1 ) {
                // This plane is thus far unique.
                unique.push_back( plane_pins );
                // Create a Plane instance for this collection of Lattices
                std::vector<const Lattice*> lattices;
                for (int ilat=0; ilat<core_.nx()*core_.ny(); ilat++) {
                    const Lattice* lat = &core_.at(ilat)[iz];
                    lattices.push_back(lat);
                }
                planes_.emplace_back(lattices, core_.nx(), core_.ny());
                unique_plane_.push_back( planes_.size() - 1 );
                first_unique_.push_back( iz );
            } else {
                // We did find a match to a previous plane. Push that ID
                unique_plane_.push_back( match_plane );
            }
            plane_pins.clear();
        } // Unique plane search

        // Put together the list of pin boundaries. For now we are treating them
        // as independent of axial plane
        x_vec_.push_back(0.0);
        real_t h_prev = 0.0;
        for( int ilatx=0; ilatx < core_.nx(); ilatx++ ) {
            const Assembly* asy = &core_.at(ilatx, 0);
            const Lattice* lat = &((*asy)[0]);
            for( auto &h: lat->hx_vec() ) {

                dx_vec_.push_back(h);
                x_vec_.push_back(h + h_prev);
                lines_.push_back( Line( Point2(h+h_prev, 0.0),
                                        Point2(h+h_prev, hy_) ) );
                h_prev += h;
            }
        }
        y_vec_.push_back(0.0);
        h_prev = 0.0;
        for( int ilaty=0; ilaty < core_.ny(); ilaty++ ) {
            const Assembly* asy = &core_.at(0, ilaty);
            const Lattice* lat = &((*asy)[0]);
            for( auto &h: lat->hy_vec() ) {
                dy_vec_.push_back(h);
                y_vec_.push_back(h + h_prev);
                lines_.push_back( Line( Point2(0.0, h+h_prev),
                                        Point2(hx_, h+h_prev) ) );
                h_prev += h;
            }
        }

        dz_vec_ = core_.dz();

        // Coarse mesh volumes.
        /**
         * \todo Really need to clean up \ref Mesh and \ref CoreMesh
         * construction. The data members of the base \ref Mesh class should
         * really be initialized by the Mesh class itself. Its hard to know the
         * necessary data to pass to the Mesh constructor before parsing the
         * core first. Going to have to come up with something.
         */

        vol_ = VecF( this->n_pin() );
        for( size_t i=0; i<this->n_pin(); i++ ) {
            auto pos = this->coarse_position( i );
            vol_[i] = dx_vec_[pos.x] * dy_vec_[pos.y] * dz_vec_[pos.z];
        }

        // Add up the number of regions and XS regions in the entire problem
        // geometry
        n_reg_   = 0;
        n_xsreg_ = 0;
        for ( auto &a: core_.assemblies() ) {
            n_reg_ += a->n_reg();
            n_xsreg_ += a->n_xsreg();
        }


        // calculate surface indices
        this->prepare_surfaces();

        return;
    } // constructor

    CoreMesh::~CoreMesh() {
        return;
    }


    const PinMeshTuple CoreMesh::get_pinmesh( Point2 &p, size_t iz,
            int &first_reg ) const {
        assert( (iz >= 0) & (iz<planes_.size()) );

        // Locate the Position of the pin
        unsigned int ix = std::lower_bound(x_vec_.begin(), x_vec_.end(), p.x ) -
            x_vec_.begin() - 1;
        unsigned int iy = std::lower_bound(y_vec_.begin(), y_vec_.end(), p.y ) -
            y_vec_.begin() - 1;

        Position pos(ix, iy, iz);

        return PinMeshTuple( pos, planes_[iz].get_pinmesh(p, first_reg) );
    }

    Position CoreMesh::pin_position( size_t ipin ) const {
        Position pos = planes_[0].pin_position( ipin % (nx_*ny_) );
        pos.z = ipin/(nx_ * ny_);
        return pos;
    }

    std::ostream& operator<<( std::ostream &os, const CoreMesh &mesh ) {
        os << "Boundary conditions: " << std::endl;
        for( int ib=0; ib<6; ib++ ) {
            os << (Surface)ib << ":\t" << mesh.bc_[ib] << std::endl;
        }
        os << std::endl;

        os << "Mesh X Pitches:" << std::endl;
        for( auto v: mesh.dx_vec_ ) {
            os << v << std::endl;
        }
        os << std::endl;

        os << "Mesh Y Pitches:" << std::endl;
        for( auto v: mesh.dy_vec_ ) {
            os << v << std::endl;
        }
        os << std::endl;

        os << "Mesh Z Pitches:" << std::endl;
        for( auto v: mesh.dz_vec_ ) {
            os << v << std::endl;
        }
        os << std::endl;

        os << "Pin Meshes: " << std::endl;
        for( const auto& pm: mesh.pin_meshes_ ) {
            os << "Mesh ID: " << pm.first << endl;
            os << *(pm.second) << std::endl;
            os << std::endl;
        }

        return os;
    }

}
