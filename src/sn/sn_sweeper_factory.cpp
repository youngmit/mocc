#include "sn_sweeper_factory.hpp"

#include <iostream>
#include <string>

#include "mocc-core/error.hpp"
#include "mocc-core/transport_sweeper.hpp"

#include "sn/sn_sweeper.hpp"
#include "sn/sn_sweeper_cdd.hpp"
#include "sn/sn_sweeper_dd.hpp"

namespace mocc {
    UP_Sweeper_t SnSweeperFactory( const pugi::xml_node &input,
            const CoreMesh &mesh ) {

        std::string equation = "dd";
        if( input.attribute("equation") ) {
            equation = input.attribute("equation").value();
        }
        std::cout << "Generating sn sweeper with equation: "
            << equation << std::endl;

        if( equation == "dd") {
            return UP_Sweeper_t( new sn::SnSweeper<sn::CellWorker_DD>( input, 
                        mesh ) );
        } else if( equation == "cdd" ) {
            return UP_Sweeper_t(
                    new sn::SnSweeper<sn::CellWorker_CDD_DD>( input,
                        mesh ) );
        } else {
            throw EXCEPT("Unrecognized equation for Sn sweeper.");
        }

        // Shouldnt ever get here. This is just to suppress warnings.
        return UP_Sweeper_t( nullptr );
    }
}
