/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "input_proc.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/file_scrubber.hpp"
#include "util/files.hpp"
#include "util/omp_guard.h"
#include "util/timers.hpp"
#include "core/angular_quadrature.hpp"
#include "core/globals.hpp"
#include "core/material_lib.hpp"
#include "core/pin_mesh.hpp"
#include "auxiliary/geometry_output.hpp"

using std::shared_ptr;

namespace {
std::string strip_extension(std::string input)
{
    std::string ret = input;
    size_t pos      = ret.rfind(".");
    return ret.substr(0, pos);
}

void apply_amendment(pugi::xml_node &node, std::string path,
                     const std::string &value);
}

namespace mocc {
InputProcessor::InputProcessor(std::vector<std::string> args)
    : timer_(RootTimer.new_timer("Input Processor", true)),
      core_mesh_(nullptr),
      solver_(nullptr),
      args_(args)
{
    std::vector<std::string> replacements;
    std::string filename      = "";
    bool good_cmd             = true;
    std::string command_error = "";
    for (size_t iarg = 1; iarg < args_.size(); iarg++) {
        std::string arg = args_[iarg];
        if (arg == "-a") {
            // Make sure that there is a next argument
            if (iarg == args_.size() - 1) {
                good_cmd = false;
            } else {
                // Make sure that the next argument isnt another flag
                if (args_[iarg + 1][0] == '-') {
                    good_cmd = false;
                }
            }

            if (!good_cmd) {
                command_error = "-a option specified without argument";
                break;
            }

            // Read the replacement string. Pre-increment is intended
            replacements.push_back(args_[++iarg]);
        } else {
            // This should be the filename
            if (filename == "") {
                filename = args_[iarg];
            } else {
                // Filename already defined
                good_cmd = false;
                std::stringstream error;
                error << "Filename appears to be multiply-defined: "
                      << filename;
                command_error = error.str();
            }
        }
    }

    // Make sure we got a filename
    if (good_cmd && filename == "") {
        good_cmd      = false;
        command_error = "No filename";
    }

    if (!good_cmd) {
        std::cerr << command_error << std::endl;

        std::cout << "Usage: mocc [-a substitution/path/attribute=value] "
                     "infile"
                  << std::endl;

        exit(EXIT_FAILURE);
    }

    case_name_ = strip_extension(filename);

    Timer &xml_timer = timer_.new_timer("XML parsing");
    xml_timer.tic();

    LogScreen << "Parsing: " << filename << std::endl;

    pugi::xml_parse_result result = doc_.load_file(filename.c_str());

    // Make sure this worked
    if (result.status != pugi::status_ok) {
        std::cerr << "XML parse error: " << result.description() << std::endl;
        throw EXCEPT("Error encountered in parsing XML file.");
    }

    // Handle the replacements specified on the command line if any
    for (const auto &replacement : replacements) {
        size_t loc = replacement.find("=");
        if (loc == std::string::npos) {
            std::stringstream msg;
            msg << "Malformed replacement in command line. Proper syntax:"
                << std::endl;
            msg << "-a path/to/attribute:new_value" << std::endl;
            throw EXCEPT(msg.str());
        }
        size_t rloc = replacement.rfind("=");
        if (rloc != loc) {
            throw EXCEPT("Malformed replacement in command line");
        }

        std::string path(replacement, 0, loc);
        std::string value(replacement, loc + 1);
        apply_amendment(doc_, path, value);
    }

    // Get problem-global stuff
    if (!doc_.child("case_name").empty()) {
        case_name_ = doc_.child("case_name").child_value();
        if (case_name_.empty()) {
            throw EXCEPT("<case_name> was provided, yet empty");
        }
    }

    global::case_name = case_name_;

    xml_timer.toc();

    

    timer_.toc();

    return;
}

void InputProcessor::process()
{
    timer_.tic();
    Timer &mesh_timer = timer_.new_timer("Core Mesh");
    mesh_timer.tic();

    // Dump the text of the input file to the log file
    std::stringstream xmlstream;
    doc_.save(xmlstream);
    LogFile << "XML input (including command-line amendments):" << std::endl;
    LogFile << " =============================================================="
               "==============="
            << std::endl;
    LogFile << xmlstream.str();
    LogFile << " =============================================================="
               "==============="
            << std::endl
            << std::endl;

    // Read the <parallel /> tag, if present. This isnt the best place for
    // this in the long term, but since we only do OpenMP right now we can
    // live with it. If MPI becomes a thing, it would be better to have a
    // standalone parallel environment class to take care of all of this.
    if (!doc_.child("parallel").empty()) {
        int n_thread =
            doc_.child("parallel").attribute("num_threads").as_int(0);

        if (n_thread < 1) {
            throw EXCEPT("Less than one thread specified in <parallel> "
                         "tag");
        }

        if (n_thread > omp_get_num_procs()) {
            Warn("More threads specified than physical "
                 "threads on this machine in <parallel> tag");
        }

        // Okay. Looking good. Tell OpenMP whats up
        omp_set_num_threads(n_thread);
    }

    // Generate the core mesh
    core_mesh_ = std::make_shared<CoreMesh>(doc_);

    mesh_timer.toc();

    Timer &solver_timer = timer_.new_timer("Solver");
    solver_timer.tic();

    // Generate a top-level solver
    solver_ = SolverFactory(doc_.child("solver"), *core_mesh_.get());

    // Perform geometry output if necessary
    if (!doc_.child("geometry_output").empty()) {
        aux::output_geometry(doc_.child("geometry_output"), *core_mesh_.get());
    }

    LogFile << std::endl;

    solver_timer.toc();
    timer_.toc();

    return;
}
};

namespace {
using namespace mocc;
void apply_amendment(pugi::xml_node &node, std::string path,
                     const std::string &value)
{
    size_t pos = path.find("/");
    if (pos == std::string::npos) {
        // No slashes. This should be either the name of the attribute in the
        // current node to alter, or the name of the child node conatining the
        // text to alter.
        //
        // Require that this be unambiguous, throwing an error if there is an
        // attribute and a child node of the same node (should probably never
        // happen anyways)
        bool has_attribute = !node.attribute(path.c_str()).empty();
        bool has_child     = !node.child(path.c_str()).empty();
        if (!(has_attribute ^ has_child)) {
            std::cout << path << std::endl;
            std::cout << has_attribute << " " << has_child << std::endl;
            throw EXCEPT("Could not find the requested attribute or child to "
                         "modify, or could not disambiguate.")
        }

        if (has_attribute) {
            if (!node.attribute(path.c_str()).set_value(value.c_str())) {
                throw EXCEPT("Failed to modify the requested attribute.");
            }
        }

        if(has_child) {
            pugi::xml_text text = node.child(path.c_str()).text();
            if(!text) {
                throw EXCEPT("No text data found at requested location.");
            }
            text = value.c_str();
        }
    } else {
        // There is a slash in there. Dig deeper
        std::string new_node_name(path, 0, pos);
        if (node.child(new_node_name.c_str()).empty()) {
            throw EXCEPT("Could not find node in command-line replacement");
        }

        std::string new_path(path, pos + 1);

        pugi::xml_node new_node = node.child(new_node_name.c_str());
        apply_amendment(new_node, new_path, value);
    }

    return;
}
}
