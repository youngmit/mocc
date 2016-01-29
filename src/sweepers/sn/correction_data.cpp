#include "correction_data.hpp"

#include <iomanip>

using std::setfill;
using std::setw;

namespace mocc {
    void CorrectionData::from_data( const pugi::xml_node &input ) {
        // First, validate the data tags. Make sure that they are the right
        // size and have cover all planes in the mesh.
        {
            int last_plane = -1;
            for( auto data = input.child("data"); data;
                data = data.next_sibling("data") )
            {
                int top_plane = data.attribute("top_plane").as_int(-1);
                if( top_plane < 0 ) {
                    throw EXCEPT("Invalid top_plane in <data />");
                }
                if( top_plane <= last_plane ) {
                    throw EXCEPT("Out-of-order or duplicate top_plane in "
                            "<data> tags" );
                }
                
                if(data.attribute("file").empty()) {
                    throw EXCEPT("No file specified.");
                }
            }
            if( last_plane != (nz_-1) ) {
                throw EXCEPT("Data do not span entire mesh.");
            }
        }

        for( int ig=0; ig<ngroup_; ig++ ) {
            for( int iang=0; iang<nang_; iang++ ) {
                
            }
        }

        return;
    }
    
    void CorrectionData::output( H5Node &file ) const {
        VecI dims;
        dims.push_back(nz_);
        dims.push_back(ny_);
        dims.push_back(nx_);

        file.create_group("/alpha_x");
        file.create_group("/alpha_y");
        file.create_group("/beta");

        for( int g=0; g<ngroup_; g++ ) {
            for( int a=0; a<nang_; a++ ) {
                {
                    std::stringstream setname;
                    setname << "/beta/" << g << "_" << a;
                 //   HDF::Write( file, setname.str(), beta, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_x/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                //    HDF::Write( file, setname.str(), alpha_x, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_y/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                   // HDF::Write( file, setname.str(), alpha_y, dims );
                }
            }
        }
        return;
    }
}
