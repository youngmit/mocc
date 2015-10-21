#pragma once

#include <map>
#include <vector>
#include <memory>
#include <cassert>

#include "pugixml.hpp"

#include "lattice.hpp"
#include "global_config.hpp"

namespace mocc {
    class Assembly {
    public:
        Assembly( const pugi::xml_node &input, 
                  const std::map<int, Lattice> &lattices );

        ~Assembly();

        unsigned int id() const {
            return id_;
        }

        /*
         * Return the number of pins along the x dimension
         */
        unsigned int nx() const {
            return lattices_[0]->nx();
        }

        /**
         * Return the number of pins along the y dimension
         */
        unsigned int ny() const {
            return lattices_[0]->ny();
        }

        /**
         * Return the number of planes
         */
        unsigned int nz() const {
            return nz_;
        }

        /**
         * Return the total height of the indexed plane in the \ref Assembly
         */
        real_t hz( unsigned int iz ) const {
            return hz_[iz];
        }

        /**
         * Return a const reference to the vector of plane heights.
         */
        const VecF& hz() const {
            return hz_;
        }
        
        // Return the total size of the assembly in the x dimension
        real_t hx() const {
            return hx_;
        }

        // Return the total size of the assembly in the y dimension
        real_t hy() const {
            return hy_;
        }

        // Return the total number of FSRs in the assembly
        unsigned int n_reg() const {
            return n_reg_;
        }

        // Return the total number of XS regions in the assembly
        unsigned int n_xsreg() const {
            return n_xsreg_;
        }

        /**
         * Return a reference to the indexed \ref Lattice.
         */
        const Lattice& operator[](unsigned int iz) const {
            assert( (iz >= 0) & (iz < lattices_.size()) );
            return *lattices_[iz];
        }

    private:
        unsigned int id_;
        unsigned int nz_;
        VecF hz_;

        real_t hx_;
        real_t hy_;

        unsigned int n_reg_;
        unsigned int n_xsreg_;

        std::vector<const Lattice*> lattices_;
    };

    typedef std::unique_ptr<Assembly> UP_Assembly_t;
}
