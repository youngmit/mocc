#include "sn_sweeper_factory_cdd.hpp"

#include "string_utils.hpp"

namespace mocc { namespace cmdo {
    CDDPair_t SnSweeperFactory_CDD( const pugi::xml_node &input,
                                    const CoreMesh &mesh)
    {
        std::string equation = "cdd";
        if( !input.attribute("equation").empty() ) {
            equation = input.attribute("equation").value();
            sanitize(equation);
        }

        if( equation != "cdd" ) {
            throw EXCEPT("Input doesn't seem to want CDD. How did you get "
                    "here?");
        }

        // Make the correction data. For now, just call from_data always. If
        // there is no data, it will return and not do anything.
        

        // Determine the type of axial treatment and create the right type of
        // sweeper.
        std::unique_ptr<SnSweeper> sweeper;
        std::shared_ptr<CorrectionData> corrections;

        std::string axial = "dd";
        if( !input.attribute("axial").empty() ) {
            axial = input.attribute("axial").value();
            sanitize(axial);
        }
        if( axial == "dd" ) {
            LogScreen << "Diamond Difference axial treatment" << std::endl;
            cmdo::SnSweeper_CDD<CellWorker_CDD_DD> *swp =
                new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_DD>( input,
                                                                  mesh );
            corrections = std::make_shared<CorrectionData>(
                mesh, swp->ang_quad().ndir()/2,
                swp->n_group());
            corrections->from_data( input );

            swp->set_corrections( corrections );
            sweeper.reset( swp );
        } else if( axial == "sc" ) {
            LogScreen << "Step Characteristics axial treatment"
                      << std::endl;
            cmdo::SnSweeper_CDD<CellWorker_CDD_SC> *swp =
                new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_SC>( input,
                                                                  mesh );
            corrections = std::make_shared<CorrectionData>(
                mesh, swp->ang_quad().ndir()/2,
                swp->n_group());
            corrections->from_data( input );

            swp->set_corrections( corrections );
            sweeper.reset( swp );
        } else if( axial == "fw" ) {
            LogScreen << "Forward Difference axial treatment" << std::endl;
            cmdo::SnSweeper_CDD<CellWorker_CDD_FW> *swp =
                new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_FW>( input,
                                                                  mesh );
            
            corrections = std::make_shared<CorrectionData>(
                mesh, swp->ang_quad().ndir()/2,
                swp->n_group());
            corrections->from_data( input );

            swp->set_corrections( corrections );
            sweeper.reset( swp );
        } else {
            throw EXCEPT("Unrecognized axial treatment in CDD.");
        }
        return CDDPair_t( std::move(sweeper), corrections );
    }
} }
