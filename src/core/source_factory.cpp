#include "source_factory.hpp"

#include <string>

#include "source_isotropic.hpp"
#include "string_utils.hpp"

enum  class ScatteringTreatment {
    P0,
    PN,
    INVALID
};

namespace mocc {
    UP_Source_t SourceFactory( const pugi::xml_node &input, 
        int n_reg, const XSMesh *xs_mesh, const ArrayB2 &flux )
    {
        if( input.empty() ) {
            throw EXCEPT("No input specified for <source>.");
        }

        UP_Source_t source;
    
        // Check scattering treatment
        ScatteringTreatment scat = ScatteringTreatment::INVALID;
        {
            std::string in_str = input.attribute("scattering").value();
            if( in_str.empty() ) {
                throw EXCEPT("No scattering treatment specified in <source />");
            }
            sanitize(in_str);
            
            if( in_str == "p0" ) {
                scat = ScatteringTreatment::P0;
            } else if ( in_str == "pn" ) {
                scat = ScatteringTreatment::PN;
            }
        }

        // Allocate/construct the source
        switch( scat ) {
        case ScatteringTreatment::P0:
            source.reset( new SourceIsotropic( n_reg, xs_mesh, flux ) );
            break;
        case ScatteringTreatment::PN:
            throw EXCEPT("Pn scattering not supported yet.");
            break;
        default:
            throw EXCEPT("Unrecognized scattering treatment in <source />");
        }

        // Apply an external source if its specified
        source->add_external( input );



        return source;
    }
} // namespace mocc
