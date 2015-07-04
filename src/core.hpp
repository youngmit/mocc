#pragma once

#include <map>
#include <vector>
#include <cassert>

#include "pugixml.hpp"

#include "arrays.hpp"
#include "assembly.hpp"

namespace mocc {
    class Core {
    public:
        Core();
        Core( const pugi::xml_node &input, 
              const std::map<int, UP_Assembly_t> &assemblies);
        ~Core();

        Assembly* at( int i ) {
            assert( 0 <= i & i < assemblies_.size() );
            return assemblies_[i];
        }

        Assembly* at( int x, int y ) const {
            assert( 0 <= x & x < nx_ );
            assert( 0 <= y & y < ny_ );
            return assemblies_[y*nx_ + x];
        }

        int nx() const {
            return nx_;
        }

        int ny() const {
            return ny_;
        }

        int nz() const {
            return assemblies_[0]->nz();
        }

    private:
        // Core dimensions (in assemblies)
        int nx_;
        int ny_;
        
        // Core dimensions (in pins)
        int npinx_;
        int npiny_;

        // 2D array of assemblies
        std::vector<Assembly*> assemblies_;
    };
}
