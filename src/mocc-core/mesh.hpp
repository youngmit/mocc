#pragma once

#include <iostream>
#include <memory>

namespace mocc {
    class Mesh {
    public:
        Mesh() { };
        Mesh( unsigned int n_reg, unsigned int n_xsreg, unsigned int nx, 
                unsigned int ny, unsigned int nz ):
            n_reg_( n_reg ),
            n_xsreg_( n_xsreg ),
            nx_( nx ),
            ny_( ny ),
            nz_( nz )
        {
std::cout << "building base mesh type" << std::endl;
            return;
        }

        unsigned int n_reg() const {
            return n_reg_;
        }
        
        unsigned int nx() const {
            return nx_;
        }

        unsigned int ny() const {
            return ny_;
        }

        unsigned int nz() const {
            return nz_;
        }

        unsigned int n_pin() const {
            return nx_*ny_*nz_;
        }

    protected:
        // Total number of FSRs in the entire geometry
        unsigned int n_reg_;
        // Total number of XS regions in the entire geometry
        unsigned int n_xsreg_;
        // Numbers of pins/planes in each dimension
        unsigned int nx_;
        unsigned int ny_;
        unsigned int nz_;
    };

    typedef std::shared_ptr<Mesh> SP_Mesh_t;
}
