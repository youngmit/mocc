#include "input_proc.hpp"

#include <memory>
#include <string>
#include <iostream>

#include "error.hpp"
#include "files.hpp"
#include "pugixml.hpp"
#include "pin_mesh.hpp"
#include "file_scrubber.hpp"
#include "material_lib.hpp"
#include "angular_quadrature.hpp"


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
            Error("Failed to open a meaningful input file. Are you sure it exists?");
        }
        
        // Generate the core mesh
        core_mesh_ = std::make_shared<CoreMesh>( doc );

        // Generate a top-level solver
        solver_ = SolverFactory( doc.child("solver"), *core_mesh_.get() );

        LogFile << endl;
        return;
    }

    


};
