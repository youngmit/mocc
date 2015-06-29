#include <iostream>
#include <exception>
#include "global_config.hpp"
#include "files.hpp"
#include "input_proc.hpp"
#include "error.hpp"


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

    // For now, assume that we are doing an eigenvalue solve. Might support FSS
    // later.

	
	// Parse the input file
	InputProc inProc(argv[1]);
}
