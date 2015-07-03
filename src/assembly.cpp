#include "assembly.hpp"

#include <string>
#include <sstream>
#include <iostream>

#include "error.hpp"
#include "string_utils.hpp"

using std::string;
using std::stringstream;

namespace mocc {
    Assembly::Assembly( const pugi::xml_node &input, 
                        const std::map<int, Lattice> &lattices ) {
        // Parse assembly ID
        id_ = input.attribute("id").as_int(-1);
        if (id_ == -1) {
            Error("Invalid assembly ID.");
        }
        
        // Parse number of planes
        nz_ = input.attribute("np").as_int(-1);
        if (nz_ == -1) {
            Error("Invalid number of planes (nz) when parsing assembly.");
        }
        
        // Parse plane heights (scalar form)
        bool scalar_hz = false;
        float hz = input.attribute("hz").as_float(0.0f);
        if (hz > 0.0f) {
            scalar_hz = true;
            // Fill the hz vector with all entries the same.
            hz_ = VecF(nz_, hz);
        }

        // Parse plane heights (array form)
        if ( auto hz_in = input.child("hz") ){
            if ( scalar_hz ) {
                // hz is over-defined
                Error("Plane heights are over-specified for assembly.");
            }
        }

        // Parse lattice IDs
        {
            string lat_str = input.child("lattices").child_value();
            stringstream inBuf( trim(lat_str) );
            while (!inBuf.eof()) {
                int lat_id;
                inBuf >> lat_id;
                if ( lattices.count(lat_id) > 0 ) {
                    lattices_.push_back( &(lattices.at(lat_id)) );
                } else {
                    Error("Unrecognized lattice ID in assembly.");
                }
            }
        }
        
        // Store lattice dimensions
        hx_ = lattices_[0]->hx();
        hy_ = lattices_[0]->hy();

        return;
    }

    Assembly::~Assembly(){
        return;
    }
}
