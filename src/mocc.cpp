#include <iostream>
#include <exception>
#include "globalconfig.hpp"
#include "files.hpp"
#include "inputproc.hpp"
#include "error.hpp"


using std::cout;
using std::cin;
using std::endl;

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
}