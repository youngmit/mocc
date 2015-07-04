#pragma once

#include <map>
#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "lattice.hpp"
#include "global_config.hpp"

namespace mocc {
    class Assembly {
    public:
        Assembly( const pugi::xml_node &input, 
                  const std::map<int, Lattice> &lattices );

        ~Assembly();

        int id() const {
            return id_;
        }

        int nx() const {
            return lattices_[0]->nx();
        }

        int ny() const {
            return lattices_[0]->ny();
        }

        int nz() const {
            return nz_;
        }

        float_t hz( int iz ) const {
            return hz_[iz];
        }
        
        float_t hx() const {
            return hx_;
        }

        float_t hy() const {
            return hy_;
        }

        const Lattice& operator[](int iz) const {
            return *lattices_[iz];
        }

    private:
        int id_;
        int nz_;
        VecF hz_;

        float_t hx_;
        float_t hy_;

        std::vector<const Lattice*> lattices_;
    };

    typedef std::unique_ptr<Assembly> UP_Assembly_t;
}
