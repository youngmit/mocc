#include "sn_sweeper_factory.hpp"

#include <iostream>
#include <memory>
#include <string>

#include "core/error.hpp"
#include "core/transport_sweeper.hpp"

#include "sn/correction_data.hpp"
#include "sn/sn_sweeper.hpp"
#include "sn/sn_sweeper_variant.hpp"
#include "sn/sn_sweeper_dd.hpp"

#include "cmdo/sn_sweeper_cdd.hpp"

using mocc::sn::UP_SnSweeper_t;

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
            return UP_SnSweeper_t( new sn::SnSweeperVariant<sn::CellWorker_DD>(
                        input, mesh ) );
        } else if( equation == "cdd" ) {
            // For now, we are assuming if we are creating a CDD sweeper from
            // this factory that it is not getting correction data from another
            // coupled sweeper. Therefore, we should make some correction
            // factors for it here.
            cmdo::SnSweeper_CDD *swp = new cmdo::SnSweeper_CDD( input, mesh );
            auto corrections = std::shared_ptr<CorrectionData>(
                new CorrectionData( mesh, swp->ang_quad().ndir()/2,
                swp->n_group()) );
            corrections->from_data( input );

            swp->set_corrections( corrections );

            return UP_SnSweeper_t( std::move(swp) );
        } else {
            throw EXCEPT("Unrecognized equation for Sn sweeper.");
        }

        // Shouldnt ever get here. This is just to suppress warnings.
        return UP_SnSweeper_t( nullptr );
    }
}
