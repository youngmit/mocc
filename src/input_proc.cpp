#include "input_proc.hpp"

#include <memory>
#include <string>
#include <iostream>

#include "pugixml.hpp"

#include "mocc-core/angular_quadrature.hpp"
#include "mocc-core/error.hpp"
#include "mocc-core/file_scrubber.hpp"
#include "mocc-core/files.hpp"
#include "mocc-core/material_lib.hpp"
#include "mocc-core/pin_mesh.hpp"

#include "auxiliary/geometry_output.hpp"


using std::cout;
using std::endl;
using std::shared_ptr;

namespace mocc{
    InputProc::InputProc(const char* filename){
        LogFile << "Processing input" << endl;
        LogFile << "Parsing: " << filename << endl;

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file( filename );

        // Make sure this worked
        if( result.status != pugi::status_ok ) {
            std::cout << "XML parse error: " << result.description()
                << std::endl;
            throw EXCEPT("Error encountered in parsing XML file.");
        }

        // Generate the core mesh
        core_mesh_ = std::make_shared<CoreMesh>( doc );

        // Generate a top-level solver
        solver_ = SolverFactory( doc.child("solver"), *core_mesh_.get() );

        // Perform geometry output if necessary
        if( !doc.child("geometry_output").empty() ) {
            aux::output_geometry( doc.child("geometry_output"),
                    *core_mesh_.get() );
        }

        LogFile << endl;
        return;
    }
};
