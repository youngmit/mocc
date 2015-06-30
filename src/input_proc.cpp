#include "input_proc.hpp"
#include "files.hpp"
#include "pugixml.hpp"
#include "pin_mesh.hpp"
#include "file_scrubber.hpp"
#include "material_lib.hpp"

#include <memory>
#include <string>
#include <iostream>

using std::cout;
using std::endl;
using std::shared_ptr;

namespace mocc{
    InputProc::InputProc(const char* filename){
        LogFile << "Processing input" << endl;
        LogFile << "Parsing: " << filename << endl;
        
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file( filename );
        
        core_mesh_ = std::make_shared<CoreMesh>( doc );
        
        LogFile << endl;
        return;
    }

    


};
