#include "pin_mesh_cyl.hpp"
#include "files.hpp"
#include "global_config.hpp"
#include "error.hpp"
#include "constants.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

using std::string;
using std::stringstream;

namespace mocc {
    PinMesh_Cyl::PinMesh_Cyl(const pugi::xml_node &input):
        PinMesh( input ) {
    	// Extract the radii and check for sanity
    	{
    		stringstream radiiIn(input.child("radii").child_value());
    		while (!radiiIn.eof()) {
    			mocc::float_t rad;
    			radiiIn >> rad;
    			xs_radii_.push_back(rad);
    			
    			if(radiiIn.fail()){
    				stringstream msg;
    				msg << "Ran into a problem reading radii for pin ID="
    				    << id_; 
    				Error(msg.str().c_str());
    			}
    		}
    		if(!radiiIn.eof()){
    			stringstream msg;
    			msg << "Dangling data in radii for pin ID="
    			    << id_;
    		}
    		// Make sure the radii are ordered
    		for (int i=0; i<xs_radii_.size()-1; i++){
    			if (xs_radii_[i] > xs_radii_[i+1]) {
    				// The radii are not ordered. Pitch a fit
    				stringstream msg;
    				msg << "Pin radii do not appear to be ordered for pin ID="
    				    << id_;
    				Error(msg.str().c_str());
    			}
    		}
    		
    		// Make sure the last radius is smaller than a half-pitch
    		if(*(xs_radii_.end()-1) > pitch_x_*0.5){
    			Error("Largest radius is too big!");
    		}
    		
    		n_xsreg_ = xs_radii_.size() + 1;
    	}
    	
    	// Read in the azimuthal subdivisions
    	{
            sub_azi_ = std::vector<int>();
    		stringstream inBuf(input.child("sub_azi").child_value());
    		int azi;
    		inBuf >> azi;
    		if(inBuf.fail()){
    			Error("Improper input to azimuthal subdivisions!");
    		}
    		sub_azi_.push_back(azi);
    		
    		// For now im only supporting one entry (same azi for all rings).
    		if(sub_azi_.size() > 1){
    			Error("Only supporting on azi type for now.");
    		}
    	}
    
    	// Read in the radial subdivisions
    	{
    		stringstream inBuf(input.child("sub_radii").child_value());
    		while (!inBuf.eof()) {
    			int sub;
    			inBuf >> sub;
    			sub_rad_.push_back(sub);
    			
    			if(inBuf.fail()){
    				stringstream msg;
    				msg << "Ran into a problem reading radial subdivisions for pin ID="
    				    << id_; 
    				Error(msg.str().c_str());
    			}
    		}
    		if(!inBuf.eof()){
    			stringstream msg;
    			msg << "Dangling data in radial subdivisions for pin ID="
    			    << id_;
    		}
    		
    		// Make sure we have the same number of radial subdivs as rings
    		if(sub_rad_.size() != n_xsreg_ - 1){
    			Error("Wrong number of radial subdivisions specified.");
    		}
    	}
    	
    	//
    	// We should be done extracting information from the XML at this point
    	//

    	// Calculate actual mesh radii. They should have equal volume within each
    	// XS ring.
    	double rxsi = 0.0;
    	double ri = 0.0;
    	for(int ixs=0; ixs<n_xsreg_ - 1; ixs++){
    		double vn = (xs_radii_[ixs]*xs_radii_[ixs] - rxsi*rxsi) / 
    		            sub_rad_[ixs];
    		
    		for(int ir=0; ir<sub_rad_[ixs]; ir++){
    			double r = sqrt(vn + ri*ri);
    			radii_.push_back(r);
    			ri = r;
    		}
    		rxsi = xs_radii_[ixs];
    	}
    	
    	n_reg_= radii_.size() * sub_azi_[0];

    	return;
    }

    PinMesh_Cyl::~PinMesh_Cyl() {
    }
}
