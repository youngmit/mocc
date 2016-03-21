#include "driver.hpp"
#include "error.hpp"

int main( int argc, char* argv[] ) {
    // Make sure we have an input file
    if(argc < 2){
        mocc::Error("No input file specified!");
    }
    return run( argv[1] );
}
