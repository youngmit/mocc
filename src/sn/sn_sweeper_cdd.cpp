#include "sn_sweeper_cdd.hpp"

#include "error.hpp"
#include "files.hpp"
#include "sn_current_worker.hpp"


using std::cout;
using std::cin;
using std::endl;

namespace mocc {
    SnSweeper_CDD::SnSweeper_CDD( const pugi::xml_node &input, 
            const CoreMesh &mesh):
        SnSweeper( input, mesh ),
        cell_worker_( mesh, ang_quad_ ),
        corrections_( nullptr )
    {
        LogFile << "Constructing a CDD Sn sweeper" << std::endl;

        if( !input.child("data").empty() ) {
            LogFile << "Located auxiliary data specification." << std::endl;
            if( !input.child("data").attribute("type").empty() ) {
                std::string data_type = 
                    input.child("data").attribute("type").value();
                if( data_type == "default" ) {
                    LogFile << "Generating default values for correction "
                        "factors." << std::endl;
                    my_corrections_.reset( new CorrectionData( n_reg_, 
                                ang_quad_.ndir(), n_group_) );
                    this->set_corrections( my_corrections_.get() );
                } else {
                    throw EXCEPT( "Unrecognized data type specified for Sn CDD "
                            "sweeper." );
                }
            } else {
                throw EXCEPT( "The <data> tag for an Sn sweeper must have a "
                        "type attrubute." );
            }
        }
    }

    void SnSweeper_CDD::sweep( int group ) {
        // Make sure we have correction factors
        if( !corrections_ ) {
            throw EXCEPT( "CDD sweeper doesn't have any correction data. Try "
                    "adding <data type=\"default\"/> in the input file." );
        }
        cell_worker_.set_group( group );

        // Store the transport cross section somewhere useful
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_[ireg] = xstr;
            }
        }
        
        flux_1g_ = flux_[std::slice(group*n_reg_, n_reg_, 1)];

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // Set the source (add upscatter and divide by 4PI)
            source_->self_scatter( group, flux_1g_, q_ );

            if( inner == n_inner_-1 && coarse_data_ ) {
                // Wipe out the existing currents
                coarse_data_->current.col( group ) = 0.0;
                this->sweep_1g<sn::Current, CellWorker_CDD>( group, 
                        cell_worker_ );
            } else {
                this->sweep_1g<sn::NoCurrent, CellWorker_CDD>( group, 
                        cell_worker_ );
            }
        }
        flux_[std::slice(group*n_reg_, n_reg_, 1)] = flux_1g_;

        return;
    }
}
