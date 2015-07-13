#pragma once

#include <map>
#include <vector>
#include <cassert>

#include "pugixml.hpp"

#include "arrays.hpp"
#include "assembly.hpp"
#include "geom.hpp"

namespace mocc {
    class Core {
    public:
        Core();
        Core( const pugi::xml_node &input, 
              const std::map<int, UP_Assembly_t> &assemblies);
        ~Core();

        const Assembly& at( unsigned int i ) const {
            assert( (0 <= i) & (i < assemblies_.size()) );
            return *assemblies_[i];
        }

        const Assembly& at( unsigned int x, unsigned int y ) const {
            assert( (0 <= x) & (x < nx_) );
            assert( (0 <= y) & (y < ny_) );
            return *assemblies_[y*nx_ + x];
        }

        const std::vector<Assembly*>& assemblies() const {
            return assemblies_;
        }

        // Return the number of assemblies along X direction
        int nx() const {
            return nx_;
        }

        // Return the number of assemblies along Y direction
        int ny() const {
            return ny_;
        }

        // Return the number of planes in the core
        int nz() const {
            return assemblies_[0]->nz();
        }

    private:
        // Core dimensions (in assemblies)
        unsigned int nx_;
        unsigned int ny_;
        
        // Core dimensions (in pins)
        unsigned int npinx_;
        unsigned int npiny_;

        // Boundaries of the lattices in the core
        VecF hx_vec_;
        VecF hy_vec_;


        // 2D array of assemblies
        std::vector<Assembly*> assemblies_;
    };
}