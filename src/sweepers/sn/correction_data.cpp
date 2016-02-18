#include "correction_data.hpp"

#include <iomanip>

#include "core/files.hpp"
#include "core/string_utils.hpp"

using std::setfill;
using std::setw;

namespace mocc {
    void CorrectionData::from_data( const pugi::xml_node &input ) {
        LogFile << "Loading CDD data from file(s)." << std::endl;
        // First, validate the data tags. Make sure that they are the right
        // size and have cover all planes in the mesh.
        {
            std::vector<bool> plane_data(nz_, false);
            for( auto data = input.child("data"); data;
                data = data.next_sibling("data") )
            {
                int top_plane = 0;
                int bot_plane = 0;
                if( !data.attribute("top_plane").empty() ) {
                    top_plane = data.attribute("top_plane").as_int();
                }
                if( !data.attribute("bottom_plane").empty() ) {
                    bot_plane = data.attribute("bottom_plane").as_int();
                }

                if( (bot_plane < 0) || (bot_plane >= (int)mesh_->nz()) ) {
                    throw EXCEPT("Invalid bottom_plane");
                }
                if( (top_plane < 0) || (top_plane >= (int)mesh_->nz()) ) {
                    throw EXCEPT("Invalid top_plane");
                }

                if(data.attribute("file").empty()) {
                    throw EXCEPT("No file specified.");
                }
                
                for( int ip=bot_plane; ip<=top_plane; ip++ ) {
                    if( plane_data[ip] ) {
                        std::stringstream msg;
                        msg << "Plane data is over-specified. Look at plane "
                            << ip;
                        throw EXCEPT(msg.str());
                    }
                    plane_data[ip] = true;
                }
            }
            LogFile << "Correction data is being specified for the following "
                "planes:" << std::endl;
            LogFile << print_range(plane_data) << std::endl;
        }

        ArrayB1 slice(mesh_->n_cell_plane());

        for( auto data = input.child("data"); data;
            data = data.next_sibling("data") )
        {
            int top_plane = 0;
            int bot_plane = 0;
            if( !data.attribute("top_plane").empty() ) {
                top_plane = data.attribute("top_plane").as_int();
            }
            if( !data.attribute("bottom_plane").empty() ) {
                bot_plane = data.attribute("bottom_plane").as_int();
            }
            H5Node h5f( data.attribute("file").value(), H5Access::READ );

            for( int ig=0; ig<ngroup_; ig++ ) {
                for( int iang=0; iang<nang_; iang++ ) {
                    // Gobble that data. Om nom nom!
                    {
                        std::stringstream path;
                        path << "/alpha_x/" << setfill('0') << setw(3) << ig <<
                            "_" << setfill('0') << setw(3) << iang;

                        h5f.read_1d( path.str(), slice );
                        for( int ip=bot_plane; ip<=top_plane; ip++ ) {
                            int stt = mesh_->plane_cell_begin(ip);
                            int stp = mesh_->plane_cell_end(ip)-1;
                            alpha_(ig, iang, blitz::Range(stt, stp), 
                                    (int)Normal::X_NORM) = slice;
                            
                        }
                    }
                    
                    {
                        std::stringstream path;
                        path << "/alpha_y/" << setfill('0') << setw(3) << ig <<
                            "_" << setfill('0') << setw(3) << iang;
                        h5f.read_1d( path.str(), slice );
                        for( int ip=bot_plane; ip<=top_plane; ip++ ) {
                            int stt = mesh_->plane_cell_begin(ip);
                            int stp = mesh_->plane_cell_end(ip)-1;
                            alpha_(ig, iang, blitz::Range(stt, stp), 
                                    (int)Normal::Y_NORM) = slice;
                        }
                    }
                    {
                        std::stringstream path;
                        path << "/beta/" << setfill('0') << setw(3) << ig << "_"
                            << setfill('0') << setw(3) << iang;
                        h5f.read_1d( path.str(), slice );
                        for( int ip=bot_plane; ip<=top_plane; ip++ ) {
                            int stt = mesh_->plane_cell_begin(ip);
                            int stp = mesh_->plane_cell_end(ip)-1;
                            beta_(ig, iang, blitz::Range(stt, stp)) = slice;
                        }
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
