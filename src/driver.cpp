#include "driver.hpp"

#include <csignal>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <sstream>

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
std::unique_ptr<InputProc> input_proc;

// Generate output from the solver
void generate_output() {
    std::string out_name = input_proc->case_name();
    out_name.append(".h5");
    H5Node outfile( out_name, H5Access::WRITE );
    solver->output( outfile );

    LogFile << std::endl;
    LogFile << "Full input:" << std::endl;

    std::stringstream filestream;

    input_proc->document().save(filestream);
    outfile.write("input_file", filestream.str());

    outfile.write("git_sha1", std::string(g_GIT_SHA1));

    std::cout << "Output written to '" << input_proc->case_name() << "'"
              << std::endl;
}

// Print the MOCC banner. Pretty!
void print_banner();

// Signal handler for SIGINT. Calls output() and quits
void int_handler(int p) {
    std::cout << "Caught SIGINT. Bailing." << std::endl;
    generate_output();
    std::exit(EXIT_FAILURE);
}


/**
 * This does the whole shebang: open and parse the input file pointed to by \p
 * file, producing a \ref mocc::Solver and \ref mocc::CoreMesh, calling \ref 
 * mocc::Solver::solve(), then calling \ref mocc::Solver::output().
 *
 * To give some rhyme to the reason for why this isn't just \c main(): This is
 * so that integration tests can link against the driver, and call \c run(),
 * then interpret the results after, without collisions of \c main(). For the
 * actual MOCC entry point, look in \c mocc.cpp, which just calls this function.
 */
int run( std::string file ) {
    std::signal( SIGINT, int_handler );

    print_banner();

    try {
        RootTimer.tic();
        
        // Set up an input processor
        input_proc.reset( new InputProc(file) );

        // Spin up the log file. We do this after we peek at the input
        // processor for a case_name tag
        StartLogFile(input_proc->case_name());

        LogScreen << "Running case: " << input_proc->case_name() << std::endl;
        LogScreen << "Using MOCC executable built with GIT SHA1: " <<
            g_GIT_SHA1 << std::endl;
        {
            time_t t;
            time(&t);
            LogScreen << "Local time: " << ctime(&t);
        }
        LogScreen << std::endl << std::endl;

#pragma omp parallel
        {
#pragma omp master
            {
                LogFile << "Running with " << omp_get_num_threads() << " treads"
                    << std::endl;
            }
        }

        // Actually process the XML input. We waited until now to do this,
        // because we want to be able to log the progress to a file, but needed
        // a case_name from the input file to be processed.
        input_proc->process();


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

    catch(Exception e) {
        cout << "Error:" << endl;
        cout << e.what();
        return 1;
    }
    return 0;
}

void print_banner() {
    std::string space = "                         ";
    std::cout << space << "01001101010011110100001101000011" << std::endl;
    std::cout << space << " __  __   _____   _____   _____" <<std::endl;
    std::cout << space << "|  \\/  | |  _  | /  __ \\ /  __ \\" << std::endl;
    std::cout << space << "| .  . | | | | | | /  \\/ | /  \\/" << std::endl;
    std::cout << space << "| |\\/| | | | | | | |     | |    " << std::endl;
    std::cout << space << "| |  | | \\ \\_/ / | \\__/\\ | \\__/ " << std::endl;
    std::cout << space << "\\_|  |_/  \\___/   \\____/  \\____/" << std::endl;
    std::cout << space << std::endl;
    std::cout << space << "01101101011011110110001101100011 " << std::endl;
}
