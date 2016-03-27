#include "sn_sweeper_factory.hpp"

#include <iostream>
#include <memory>
#include <string>


#include "core/error.hpp"
#include "core/string_utils.hpp"
#include "core/transport_sweeper.hpp"

#include "sn/sn_sweeper.hpp"
#include "sn/sn_sweeper_variant.hpp"
#include "sn/sn_sweeper_dd.hpp"

#include "cmdo/correction_data.hpp"
#include "cmdo/sn_sweeper_factory_cdd.hpp"
#include "cmdo/sn_sweeper_cdd.hpp"

using mocc::sn::UP_SnSweeper_t;

namespace mocc {
    UP_SnSweeper_t SnSweeperFactory( const pugi::xml_node &input,
            const CoreMesh &mesh ) {
        using namespace sn;
        using namespace cmdo;

        std::string equation = "dd";
        if( input.attribute("equation") ) {
            equation = input.attribute("equation").value();
        }
        std::cout << "Generating sn sweeper with equation: "
            << equation << std::endl;

        if( equation == "dd") {
            std::string axial = "dd";
            if( !input.attribute("axial").empty() ) {
                axial = input.attribute("axial").value();
            }
            if( axial == "dd" ) {
                LogScreen << "Using Diamond Difference axial treatment"
                    << std::endl;
                return UP_SnSweeper_t
                    (
                        new sn::SnSweeperVariant<sn::CellWorker_DD>(
                            input, mesh )
                    );
            } else if( axial == "sc" ) {
                LogScreen << "Using Step Characteristics axial treatment"
                    << std::endl;
                    return UP_SnSweeper_t
                    (
                        new sn::SnSweeperVariant<sn::CellWorker_DD_SC>(
                            input, mesh )
                    );
            } else {
                throw EXCEPT("Unsupported axial treatment");
            }
        } else if( equation == "cdd" ) {
            // Defer to the CDD factory, but discard the correction data
            auto cdd_pair = SnSweeperFactory_CDD( input, mesh );
            return std::move(cdd_pair.first);
        } else {
            throw EXCEPT("Unrecognized equation for Sn sweeper.");
        }

        // Shouldnt ever get here. This is just to suppress warnings.
        return UP_SnSweeper_t( nullptr );
    }
}
