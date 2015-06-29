#include "inputproc.hpp"
#include "files.hpp"
#include "pugixml.hpp"
#include "pinmesh.hpp"
#include "filescrubber.hpp"
#include "materiallib.hpp"

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
    pugi::xml_parse_result result = doc.load_file(filename);
    
    // Parse meshes
    for (pugi::xml_node mesh = doc.child("mesh"); mesh; mesh = mesh.next_sibling("mesh")){
        LogFile << "Parsing new pin mesh: ID=" 
        << mesh.attribute("id").value() << endl;
        SP_PinMesh pm = PinMeshFactory(mesh);
        cout << "pin mesh id: " << pm->id() << endl;
        m_pinMeshes.insert(std::pair<int, SP_PinMesh>(pm->id(), pm));
    }
    
    // Parse Material Library
    std::string matLibName = 
        doc.child("material_lib").attribute("path").value();
    cout << "Found material library specification: " << matLibName << endl;
    FileScrubber matLibFile(matLibName.c_str(), "!");
    MaterialLib matLib(matLibFile);
    
    // Parse material IDs
    for(pugi::xml_node mat = doc.child("material_lib").child("material"); 
            mat; mat = mat.next_sibling("material")){
        cout << mat.attribute("id").value() << " "
             << mat.attribute("name").value() << endl;
        matLib.assignID(mat.attribute("id").as_int(),
                             mat.attribute("name").value());
    }
    
    // Parse pins
    
    LogFile << endl;
    return;
}


};
