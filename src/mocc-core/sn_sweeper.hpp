#pragma once

#include "pugixml.hpp"

#include "global_config.hpp"
#include "transport_sweeper.hpp"

namespace mocc {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh );

        ~SnSweeper() { }

        void sweep( int group ) {}

        void initialize() { }
        void get_pin_flux( int ig, VecF& flux ) const { }

        void calc_fission_source( float_t k, 
                ArrayX& fission_source) const { }

        void output( H5File& file ) const { }
    private:
    };
}
