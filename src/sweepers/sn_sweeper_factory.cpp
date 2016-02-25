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
        using namespace sn;
        using namespace cmdo;

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
            // Make the sweeper
            std::string axial = "dd";
            if( !input.attribute("axial").empty() ) {
                axial = input.attribute("axial").value();
            }
            LogScreen << "Using ";
            if( axial == "dd" ) {
                LogScreen << "Diamond Difference axial treatment" << std::endl;
                cmdo::SnSweeper_CDD<CellWorker_CDD_DD> *swp =
                    new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_DD>( input,
                                                                      mesh );
                auto corrections = std::shared_ptr<CorrectionData>(
                    new CorrectionData( mesh, swp->ang_quad().ndir()/2,
                    swp->n_group()) );
                corrections->from_data( input );
                swp->set_corrections( corrections );
                return UP_SnSweeper_t( std::move(swp) );
            } else if( axial == "sc" ) {
                LogScreen << "Step Characteristics axial treatment"
                          << std::endl;
                cmdo::SnSweeper_CDD<CellWorker_CDD_SC> *swp =
                    new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_SC>( input,
                                                                      mesh );
                auto corrections = std::shared_ptr<CorrectionData>(
                    new CorrectionData( mesh, swp->ang_quad().ndir()/2,
                    swp->n_group()) );
                corrections->from_data( input );
                swp->set_corrections( corrections );
                return UP_SnSweeper_t( std::move(swp) );

            } else if( axial == "fw" ) {
                LogScreen << "Forward Difference axial treatment" << std::endl;
                cmdo::SnSweeper_CDD<CellWorker_CDD_FW> *swp =
                    new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_FW>( input,
                                                                      mesh );
                auto corrections = std::shared_ptr<CorrectionData>(
                    new CorrectionData( mesh, swp->ang_quad().ndir()/2,
                    swp->n_group()) );
                corrections->from_data( input );
                swp->set_corrections( corrections );
                return UP_SnSweeper_t( std::move(swp) );
            } else {
                throw EXCEPT("Unrecognized axial treatment in CDD.");
            }

        } else {
            throw EXCEPT("Unrecognized equation for Sn sweeper.");
        }

        // Shouldnt ever get here. This is just to suppress warnings.
        return UP_SnSweeper_t( nullptr );
    }
}
