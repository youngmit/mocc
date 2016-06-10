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

#include "driver.hpp"

#include <csignal>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>

#include "pugixml.hpp"

#include "core/core_mesh.hpp"
#include "core/error.hpp"
#include "core/files.hpp"
#include "core/global_config.hpp"
#include "core/h5file.hpp"
#include "core/solver.hpp"
#include "core/timers.hpp"
#include "core/transport_sweeper.hpp"

#include "git_SHA1.hpp"
#include "input_proc.hpp"

using std::cout;
using std::cin;
using std::endl;

using namespace mocc;

// The global top-level solver
SP_Solver_t solver;

// Global core mesh
SP_CoreMesh_t mesh;

// Input processor
std::unique_ptr<InputProcessor> input_proc;

// Generate output from the solver
void generate_output()
{
    std::string out_name = input_proc->case_name();
    out_name.append(".h5");
    H5Node outfile(out_name, H5Access::WRITE);
    solver->output(outfile);

    LogFile << std::endl;
    LogFile << "Full input:" << std::endl;

    std::stringstream filestream;

    input_proc->document().save(filestream);
    outfile.write("input_file", filestream.str());

    outfile.write("git_sha1", std::string(g_GIT_SHA1));

    if (Warnings.size() > 0) {
        if (Warnings.size() == 1) {
            std::cout << "There was ";
        }
        else {
            std::cout << "There were ";
        }
        std::cout << Warnings.size();
        if (Warnings.size() == 1) {
            std::cout << " warning:" << std::endl;
        }
        else {
            std::cout << " warnings:" << std::endl;
        }
        for (const auto &warning : Warnings) {
            std::cout << "\t" << warning.second << std::endl;
        }
    }

    std::cout << "Output written to '" << input_proc->case_name() << "'"
              << std::endl;
}

// Print the MOCC banner. Pretty!
void print_banner();

// Signal handler for SIGINT. Calls output() and quits
void int_handler(int p)
{
    std::cout << "Caught SIGINT. Bailing." << std::endl;
    generate_output();
    std::exit(EXIT_FAILURE);
}

/**
 * This does the whole shebang: parse command line, open and parse the input
 * file, producing a \ref mocc::Solver and \ref mocc::CoreMesh, calling \ref
 * mocc::Solver::solve(), then calling \ref mocc::Solver::output().
 */
int run(int argc, char *argv[])
{
    std::vector<std::string> args;

    for (int i = 0; i < argc; i++) {
        args.push_back(argv[i]);
    }

    return run(args);
}

int run(const std::vector<std::string> &args)
{
    std::signal(SIGINT, int_handler);

    print_banner();

    try {
        RootTimer.tic();

        // Set up an input processor
        input_proc.reset(new InputProcessor(args));

        // Spin up the log file. We do this after we peek at the input
        // processor for a case_name tag
        StartLogFile(input_proc->case_name());

        LogScreen << "Running case: " << input_proc->case_name() << std::endl;
        LogScreen << "Using MOCC executable built with GIT SHA1: " << g_GIT_SHA1
                  << std::endl;
        {
            time_t t;
            time(&t);
            LogScreen << "Local time: " << ctime(&t);
        }
        LogScreen << std::endl << std::endl;

        // Actually process the XML input. We waited until now to do this,
        // because we want to be able to log the progress to a file, but needed
        // a case_name from the input file to be processed.
        input_proc->process();

#pragma omp parallel
        {
#pragma omp master
            {
                LogScreen << "Running with " << omp_get_num_threads()
                          << " treads" << std::endl;
            }
        }

        // Get an SP to the core mesh
        mesh = input_proc->core_mesh();
        LogFile << *mesh << endl;

        // Pull a shared pointer to the top-level solver and make it go
        solver = input_proc->solver();
        solver->solve();

        // Output stuff
        generate_output();

        RootTimer.toc();
        std::cout << RootTimer << std::endl;
        RootTimer.print(LogFile);

        StopLogFile();
    }

    catch (Exception e) {
        cout << "Error:" << endl;
        cout << e.what();
        return 1;
    }
    return 0;
}

void print_banner()
{
    std::string space = "                         ";
    std::cout << space << "01001101010011110100001101000011" << std::endl;
    std::cout << space << " __  __   _____   _____   _____" << std::endl;
    std::cout << space << "|  \\/  | |  _  | /  __ \\ /  __ \\" << std::endl;
    std::cout << space << "| .  . | | | | | | /  \\/ | /  \\/" << std::endl;
    std::cout << space << "| |\\/| | | | | | | |     | |    " << std::endl;
    std::cout << space << "| |  | | \\ \\_/ / | \\__/\\ | \\__/\\ "
              << std::endl;
    std::cout << space << "\\_|  |_/  \\___/   \\____/  \\____/" << std::endl;
    std::cout << space << std::endl;
    std::cout << space << "01101101011011110110001101100011 " << std::endl;
}
