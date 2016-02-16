#include "input_proc.hpp"

#include <memory>
#include <string>
#include <iostream>

#include "pugixml.hpp"

#include "core/angular_quadrature.hpp"
#include "core/error.hpp"
#include "core/file_scrubber.hpp"
#include "core/files.hpp"
#include "core/material_lib.hpp"
#include "core/pin_mesh.hpp"
#include "core/timers.hpp"

#include "auxiliary/geometry_output.hpp"


using std::cout;
using std::endl;
using std::shared_ptr;

namespace mocc{
    InputProc::InputProc(const char* filename) {
        Timer &timer = RootTimer.new_timer( "Initialization" );
        timer.tic();
        Timer &input_timer = timer.new_timer( "Input processing" );
        input_timer.tic();

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

        input_timer.toc();

        Timer &mesh_timer = timer.new_timer("Core Mesh");
        mesh_timer.tic();

        // Generate the core mesh
        core_mesh_ = std::make_shared<CoreMesh>( doc );

        mesh_timer.toc();

        // Generate a top-level solver
        solver_ = SolverFactory( doc.child("solver"), *core_mesh_.get() );

        // Perform geometry output if necessary
        if( !doc.child("geometry_output").empty() ) {
            aux::output_geometry( doc.child("geometry_output"),
                    *core_mesh_.get() );
        }

        LogFile << endl;
        timer.toc();
        return;
    }
};
