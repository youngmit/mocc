#include "core.hpp"

#include <iostream>
#include <sstream>
#include <string>

#include "error.hpp"
#include "string_utils.hpp"

namespace mocc {
    Boundary bc_parse( const pugi::xml_node& input, const char* surf ) {
        std::string in = input.attribute(surf).value();
        if ( in == "vacuum" ) {
            return Boundary::VACUUM;
        } else if ( in == "reflect" ) {
            return Boundary::REFLECT;
        } else {
            return Boundary::INVALID;
        }
    }

    Core::Core() {
        nx_ = 0;
        ny_ = 0;
    }

    Core::Core( const pugi::xml_node &input,
                const std::map<int, UP_Assembly_t> &assemblies):
        nx_( input.attribute("nx").as_int(0) ),
        ny_( input.attribute("ny").as_int(0) ),
        bc_(6, Boundary::INVALID)
    {
        // Make sure that we read the a proper ID
        if( (nx_ < 1) | (ny_ < 1) ) {
            Error("Invalid core dimensions.");
        }

        // Read in the boundary conditions
        bc_[(int)Surface::NORTH]  = bc_parse( input, "north" );
        bc_[(int)Surface::SOUTH]  = bc_parse( input, "south" );
        bc_[(int)Surface::EAST]   = bc_parse( input, "east" );
        bc_[(int)Surface::WEST]   = bc_parse( input, "west" );
        bc_[(int)Surface::TOP]    = bc_parse( input, "top" );
        bc_[(int)Surface::BOTTOM] = bc_parse( input, "bottom" );

        for( int i=0; i<6; i++ ) {
            if( bc_[i] == Boundary::INVALID ) {
                throw EXCEPT("Not all boundary conditions properly specified.");
            }
        }

        // Read in the assembly IDs
        std::string asy_str = input.child_value();
        std::stringstream inBuf( trim(asy_str) );
        std::vector<int> asy_vec;
        for ( unsigned int i=0; i<nx_*ny_; ++i ) {
            int id;
            inBuf >> id;
            asy_vec.push_back( id );
            if (inBuf.fail()) {
                Error("Trouble reading assembly IDs in core specification.");
            }
        }
        // Store references to the assemblies in a 2D array. Make sure to flip
        // the y-index to get it into lower-left origin
        assemblies_.resize( nx_*ny_ );
        int iasy=0;
        for ( unsigned int iy=0; iy<ny_; iy++ ) {
            unsigned int row = ny_ - iy - 1;
            for ( unsigned int ix=0; ix<nx_; ix++ ) {
                unsigned int col = ix;
                int asy_id = asy_vec[iasy++];
                try {
                    Assembly* asy_p =
                        assemblies.at( asy_id ).get();
                    assemblies_[row*nx_ + col] = asy_p;
                } catch(std::out_of_range) {
                    Error("Failed to locate assembly in core specification.");
                }
            }
        }

        // Check to make sure that the assemblies all fit together
        // Main things to check are that they all have the same number of planes
        // and that the heights of each plane match
        unsigned int nz = (*assemblies_.begin())->nz();
        for (auto asy = assemblies_.begin(); asy != assemblies_.end(); ++asy) {
            if ( (*asy)->nz() != nz ) {
                Error("Assemblies in the core have incompatible numbers of "
                        "planes."); }
        }

        for( unsigned int i=0; i<nz; i++ ) {
            real_t hz = (*assemblies_.begin())->hz(i);
            for( auto asy = assemblies_.begin();
                 asy != assemblies_.end(); ++asy ) {
                if ((*asy)->hz(i) != hz) {
                    Error("Assemblies have incompatible plane heights in "
                            "core.");
                }
            }
        }

        // Get the total number of pins along each dimension
        npinx_ = 0;
        for( unsigned int i=0; i<nx_; i++) {
            npinx_ += this->at(i, 0).nx();
        }
        npiny_ = 0;
        for( unsigned int i=0; i<ny_; i++) {
            npiny_ += this->at(0, i).ny();
        }

        // Store the x and y boundaries of the assemblies
        real_t prev = 0.0;
        for( unsigned int ix=0; ix<nx_; ix++) {
            hx_vec_.push_back(prev+this->at(ix, 0).hx());
            prev = hx_vec_[ix];
        }
        prev = 0.0;
        for( unsigned int iy=0; iy<ny_; iy++) {
            hy_vec_.push_back(prev+this->at(0, iy).hy());
            prev = hy_vec_[iy];
        }
    }

    Core::~Core() {
    }

    Core ParseCores( const pugi::xml_node &input, 
            const std::mat<int, UP_Assembly_t> assemblies ) {
        Core core;
        int n_core_enabled = 0;
        for( auto core_xml = input.child("core"); core_xml;
                core_xml = core_xml.next_sibling("core") ) {
            bool core_enabled = true;
            if( !core_xml.attribute("enabled").empty() ) {
                core_enabled = core_xml.attribute("enabled").as_bool();
            }
            if( core_enabled ) {
                core = Core(core_xml, assemblies_);
                n_core_enabled++;
            }
        }

        if( n_core_enabled == 0 ) {
            throw EXCEPT("No enabled core specifications.");
        }

        if( n_core_enabled > 1 ) {
            throw EXCEPT("More than one enabled core specification found. Tell "
                    "me which one to use");
        }

        return core;
    }
}
