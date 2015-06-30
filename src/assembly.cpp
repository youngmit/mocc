#include "assembly.hpp"

#include "error.hpp"

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

        
        return;
    }
}
