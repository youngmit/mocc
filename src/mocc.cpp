#include <iostream>
#include <exception>

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

int main(int argc, char* argv[]){
	cout << "Ahoy!" << endl;
	cout << "cmd: " << argv[0] << endl;
	
	// Make sure we have an input file
	if(argc < 2){
		Error("No input file specified!");
	}
	
	// Spin up the log file. For now, just use the name of the input file.
	StartLogFile(argv[1]);
	
	LogFile << "Welcome to " << PROG_NAME << "!" << std::endl << std::endl;

	// Parse the input file
	InputProc inProc(argv[1]);

    // Pull a shared pointer to the top-level solver and make it go
    SP_Solver_t solver = inProc.solver();
    solver->solve();

    // Output stuff
    const TransportSweeper* sweeper_p = solver->sweeper();

    VecF flux;
    sweeper_p->get_pin_flux(0, flux);
    
    H5File outfile("out.h5");

    outfile.write("foo", flux, VecI {6, 9});

    StopLogFile();
}
