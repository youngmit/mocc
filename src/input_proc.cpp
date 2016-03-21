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

std::string strip_extension( std::string input ) {
    std::string ret = input;
    size_t pos = ret.rfind(".");
    return ret.substr(0, pos);
}

namespace mocc{
    InputProc::InputProc(std::string filename):
        timer_(RootTimer.new_timer("Input Processor", true)),
        core_mesh_(nullptr),
        solver_(nullptr),
        case_name_(strip_extension(filename))
    {

        Timer &xml_timer = timer_.new_timer("XML parsing");
        xml_timer.tic();

        LogFile << "Processing input" << endl;
        LogFile << "Parsing: " << filename << endl;

        pugi::xml_parse_result result = doc_.load_file( filename.c_str() );

        // Make sure this worked
        if( result.status != pugi::status_ok ) {
            std::cout << "XML parse error: " << result.description()
                << std::endl;
            throw EXCEPT("Error encountered in parsing XML file.");
        }


        // Get problem-global stuff
        if( !doc_.child("case_name").empty() ) {
            case_name_ = doc_.child("case_name").child_value();
            if( case_name_.empty() ) {
                throw EXCEPT("<case_name> was provided, yet empty");
            }
        }

        xml_timer.toc();
        timer_.toc();

        return;
    }

    void InputProc::process() {
        timer_.tic();
        Timer &mesh_timer = timer_.new_timer("Core Mesh");
        mesh_timer.tic();

        // Generate the core mesh
        core_mesh_ = std::make_shared<CoreMesh>( doc_ );

        mesh_timer.toc();

        Timer &solver_timer = timer_.new_timer("Solver");
        solver_timer.tic();

        // Generate a top-level solver
        solver_ = SolverFactory( doc_.child("solver"), *core_mesh_.get() );

        // Perform geometry output if necessary
        if( !doc_.child("geometry_output").empty() ) {
            aux::output_geometry( doc_.child("geometry_output"),
                    *core_mesh_.get() );
        }

        LogFile << endl;

        solver_timer.toc();
        timer_.toc();

        return;
    }
};
