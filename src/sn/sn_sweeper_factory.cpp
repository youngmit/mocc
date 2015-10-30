#include "sn_sweeper_factory.hpp"

#include <iostream>
#include <string>

#include "error.hpp"
#include "sn_sweeper_cdd.hpp"
#include "sn_sweeper_dd.hpp"

namespace mocc {
    UP_SnSweeper_t SnSweeperFactory( const pugi::xml_node &input, 
            const CoreMesh &mesh ) {

        std::string equation = "dd";
        if( input.attribute("equation") ) {
            equation = input.attribute("equation").value();
        }
        std::cout << "Generating sn sweeper with equation: " 
            << equation << std::endl;

        if( equation == "dd") {
            return UP_SnSweeper_t( new SnSweeper_DD( input, mesh ) );
        } else if( equation == "cdd" ) {
            return UP_SnSweeper_t( new SnSweeper_CDD( input, mesh ) );
        } else {
            throw EXCEPT("Unrecognized equation for Sn sweeper.");
        }

        // Shouldnt ever get here. This is just to suppress warnings.
        return UP_SnSweeper_t( nullptr );
    }
}
