#include "core_mesh.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include "files.hpp"
#include "error.hpp"
#include "string_utils.hpp"


using std::cout;
using std::endl;
using std::stringstream;


namespace mocc {
    CoreMesh::CoreMesh(pugi::xml_node &input) {
        // Parse meshes
        for (pugi::xml_node mesh = input.child( "mesh" );
             mesh; 
             mesh = mesh.next_sibling( "mesh" )) {
            LogFile << "Parsing new pin mesh: ID=" 
            << mesh.attribute( "id" ).value() << endl;
            UP_PinMesh_t pm( PinMeshFactory( mesh ) );
            int id = pm->id();
            pin_meshes_.emplace(id, std::move( pm ) );
        }
        
        // Parse Material Library
        std::string matLibName = 
            input.child( "material_lib" ).attribute( "path" ).value();
        cout << "Found material library specification: " << matLibName << endl;
        FileScrubber matLibFile( matLibName.c_str(), "!" );
        mat_lib_ = MaterialLib( matLibFile );
        
        // Parse material IDs
        for (pugi::xml_node mat = input.child( "material_lib" ).child( "material" ); 
                mat; mat = mat.next_sibling( "material" )){
            cout << mat.attribute( "id" ).value() << " "
                 << mat.attribute( "name" ).value() << endl;
            mat_lib_.assignID( mat.attribute( "id" ).as_int(),
                               mat.attribute( "name" ).value() );
        }
        
        // Parse pins
        for ( pugi::xml_node pin = input.child( "pin" ); pin; 
                pin = pin.next_sibling( "pin" ) ) {
            
            // Construct the pin and add to the map
            UP_Pin_t pin_p( new Pin( pin, pin_meshes_ ) );
            pins_.emplace( pin_p->id(), std::move(pin_p) );
        }

        // Parse lattices
        for ( pugi::xml_node lat = input.child( "lattice" ); lat;
                lat = lat.next_sibling( "lattice" )) {
            Lattice lattice( lat, pins_ );
            lattices_.emplace( lattice.id(), lattice );
        }

        // Parse assemblies
        for ( pugi::xml_node asy = input.child("assembly"); asy;
                asy = asy.next_sibling("assembly") ) {
            int asy_id = asy.attribute("id").as_int();
            UP_Assembly_t asy_p( new Assembly( asy, lattices_ ) );
            assemblies_.emplace( asy_p->id(), std::move(asy_p) );
        }

        // Parse core
        core_ = Core( input.child("core"), assemblies_ );

        nx_ = core_.nx();
        ny_ = core_.ny();
        nz_ = core_.nz();
        nasy_ = nx_*ny_;

        // Calculate the total core dimensions
        hx_ = 0.0;
        for ( int ix=0; ix<core_.nx(); ix++ ) {
            hx_ += core_.at(ix, 0)->hx();
        }
        hy_ = 0.0;
        for ( int iy=0; iy<core_.ny(); iy++ ) {
            hy_ += core_.at(iy, 0)->hy();
        }

        // Determine the set of geometricaly-unique axial planes
        std::vector< VecI > unique;
        VecI plane_pins;
        for ( int iz=0; iz<nz_; iz++) {
            // Form a list of all pin meshes in the core plane iz
            for ( int iasy=0; iasy<nasy_; iasy++ ) {
                Assembly* asy = core_.at(iasy);
                for ( auto pin=(*asy)[iz].begin(); 
                    pin     !=(*asy)[iz].end(); ++pin ) {
                    plane_pins.push_back((*pin)->mesh_id());
                }
            }
            // Check against current list of unique planes
            int match_plane = -1;
            for( int iiz=0; iiz<unique.size(); iiz++ ) {
                for( int ip=0; ip<plane_pins.size(); ip++ ) {
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
                unique_plane_.push_back( iz );
                first_unique_.push_back( iz );
            } else {
                // We did find a match to a previous plane. Push that ID
                unique_plane_.push_back( first_unique_[match_plane] );
            }
            plane_pins.clear();
        } // Unique plane search
        return;
    }

    CoreMesh::~CoreMesh() {
        return;
    }
}
