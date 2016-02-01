#include "correction_data.hpp"

#include <iomanip>

#include "core/files.hpp"

using std::setfill;
using std::setw;

namespace mocc {
    void CorrectionData::from_data( const pugi::xml_node &input ) {
        LogFile << "Loading CDD data from file(s)." << std::endl;
        // First, validate the data tags. Make sure that they are the right
        // size and have cover all planes in the mesh.
        {
            int last_plane = -1;
            for( auto data = input.child("data"); data;
                data = data.next_sibling("data") )
            {
                int top_plane = data.attribute("top_plane").as_int(0);
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
                last_plane = top_plane;
            }
            if( last_plane != (nz_-1) ) {
                throw EXCEPT("Data do not span entire mesh.");
            }
        }

        ArrayB1 slice(nreg_);

        for( auto data = input.child("data"); data;
            data = data.next_sibling("data") )
        {
            H5Node h5f( data.attribute("file").value(), H5Access::READ );

            for( int ig=0; ig<ngroup_; ig++ ) {
                for( int iang=0; iang<nang_; iang++ ) {
                    // Gobble that data. Om nom nom!
                    {
                        std::stringstream path;
                        path << "/alpha_x/" << setfill('0') << setw(3) << ig <<
                            "_" << setfill('0') << setw(3) << iang;

                        h5f.read_1d( path.str(), slice );
                        alpha_(ig, iang, blitz::Range::all(), 
                                (int)Normal::X_NORM) = slice;
                    }
                    
                    {
                        std::stringstream path;
                        path << "/alpha_y/" << setfill('0') << setw(3) << ig <<
                            "_" << setfill('0') << setw(3) << iang;
                        h5f.read_1d( path.str(), slice );
                        alpha_(ig, iang, blitz::Range::all(), 
                                (int)Normal::Y_NORM) = slice;
                    }
                    {
                        std::stringstream path;
                        path << "/beta/" << setfill('0') << setw(3) << ig << "_"
                            << setfill('0') << setw(3) << iang;
                        h5f.read_1d( path.str(), slice );
                        beta_(ig, iang, blitz::Range::all()) = slice;
                    }
                }
            }
        }

        return;
    }
    
    void CorrectionData::output( H5Node &file ) const {
        VecI dims;
        dims.push_back(nz_);
        dims.push_back(ny_);
        dims.push_back(nx_);
        int n = nx_*ny_*nz_;

        file.create_group("/alpha_x");
        file.create_group("/alpha_y");
        file.create_group("/beta");

        // Declare slice storage
        ArrayB1 slice(n);

        for( int g=0; g<ngroup_; g++ ) {
            for( int a=0; a<nang_; a++ ) {
                {
                    slice = beta_(g, a, blitz::Range::all());
                    std::stringstream setname;
                    setname << "/beta/" << setfill('0') << setw(3) << g 
                            << "_"      << setfill('0') << setw(3) << a;
                    file.write( setname.str(), slice, dims );
                }

                {
                    slice = alpha_(g, a, blitz::Range::all(),
                            (int)Normal::X_NORM);
                    std::stringstream setname;
                    setname << "/alpha_x/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    file.write( setname.str(), slice, dims );
                }

                {
                    slice = alpha_(g, a, blitz::Range::all(),
                            (int)Normal::Y_NORM);
                    std::stringstream setname;
                    setname << "/alpha_y/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    file.write( setname.str(), slice, dims );
                }
            }
        }
        return;
    }
}
