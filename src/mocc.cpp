#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>
#include <string>

#include "global_config.hpp"
#include "files.hpp"
#include "input_proc.hpp"
#include "core_mesh.hpp"
#include "error.hpp"
#include "h5file.hpp"
#include "transport_sweeper.hpp"


using std::cout;
using std::cin;
using std::endl;

using namespace mocc;

// The global top-level solver
SP_Solver_t solver;

// Global core mesh
SP_CoreMesh_t mesh;

int main(int argc, char* argv[]){
	// Make sure we have an input file
	if(argc < 2){
		Error("No input file specified!");
	}
	
    try {
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

	    // Spin up the log file. For now, just use the name of the input file.
	    StartLogFile(argv[1]);
	    

	    // Parse the input file
	    InputProc inProc(argv[1]);

        // Get an SP to the core mesh
        mesh = inProc.core_mesh();

        // Pull a shared pointer to the top-level solver and make it go
        solver = inProc.solver();
        solver->solve();

        // Output stuff
        HDF::H5File outfile("out.h5");
        solver->output( outfile.get() );

        StopLogFile();
    }

    catch(Exception e) {
        cout << "Error:" << endl;
        cout << e.what();
        return 1;
    }
}
