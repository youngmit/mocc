#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>

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

void output();

int main(int argc, char* argv[]){
	// Make sure we have an input file
	if(argc < 2){
		Error("No input file specified!");
	}
	
	// Spin up the log file. For now, just use the name of the input file.
	StartLogFile(argv[1]);
	
	LogFile << "Welcome to " << PROG_NAME << "!" << std::endl << std::endl;

	// Parse the input file
	InputProc inProc(argv[1]);

    // Get an SP to the core mesh
    mesh = inProc.core_mesh();

    // Pull a shared pointer to the top-level solver and make it go
    solver = inProc.solver();
    solver->solve();

    // Output stuff
    output();

    StopLogFile();
}

void output(){
    const TransportSweeper* sweeper_p = solver->sweeper();

    // Get core dimensions from the mesh
    const int nx = mesh->nx();
    const int ny = mesh->ny();
    const int nz = mesh->nz();
    VecI dims;
    dims.push_back(nz);
    dims.push_back(ny);
    dims.push_back(nx);

    H5File outfile("out.h5");



    // Make a group in the file to store the flux
    outfile.mkdir("/flux");

    // Provide energy group upper bounds
    VecF eubounds = sweeper_p->eubounds();
    outfile.write("/eubounds", eubounds, VecI(1, eubounds.size()));
    outfile.write("/ng", eubounds.size());

    for( int ig=0; ig<sweeper_p->n_grp(); ig++ ) {
        VecF flux;
        sweeper_p->get_pin_flux(ig, flux);


        std::stringstream setname;
        setname << "/flux/" << std::setfill('0') << std::setw(3) << ig+1;

        outfile.write(setname.str(), flux, dims);
    }
    
}
