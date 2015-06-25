#include "pinmesh_cyl.hpp"
#include "files.hpp"
#include "globalconfig.hpp"
#include "error.hpp"
#include "constants.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

using std::string;
using std::stringstream;
using std::cout;
using std::endl;


PinMesh_Cyl::PinMesh_Cyl(const pugi::xml_node &input){
	// Extract pin id
	{
		stringstream inBuf(input.attribute("id").value());
		inBuf >> m_id;
		if(inBuf.fail()){
			Error("Failed to read pin ID.");
		}
		if(!inBuf.eof()){
			Warn("Dangling data after pin ID.");
		}
	}
	
	// Extract pitch
	{
		stringstream inBuf(input.attribute("pitch").value());
		inBuf >> m_pitchX;
		// Just treat square pitch for now
		m_pitchY = m_pitchX;
		if(inBuf.fail()){
			Error("Failed to read pin pitch.");
		}
		if(!inBuf.eof()){
			Warn("Dangling data after pin pitch.");
		}	
	}
	
	// Extract the radii and check for sanity
	{
		stringstream radiiIn(input.child("radii").child_value());
		while (!radiiIn.eof()) {
			mocc::float_t rad;
			radiiIn >> rad;
			m_radiiXS.push_back(rad);
			
			if(radiiIn.fail()){
				stringstream msg;
				msg << "Ran into a problem reading radii for pin ID="
				    << m_id; 
				Error(msg.str().c_str());
			}
		}
		if(!radiiIn.eof()){
			stringstream msg;
			msg << "Dangling data in radii for pin ID="
			    << m_id;
		}
		// Make sure the radii are ordered
		for (int i=0; i<m_radiiXS.size()-1; i++){
			if (m_radiiXS[i] > m_radiiXS[i+1]) {
				// The radii are not ordered. Pitch a fit
				stringstream msg;
				msg << "Pin radii do not appear to be ordered for pin ID="
				    << m_id;
				Error(msg.str().c_str());
			}
		}
		
		// Make sure the last radius is smaller than a half-pitch
		if(*(m_radiiXS.end()-1) > m_pitchX*0.5){
			Error("Largest radius is too big!");
		}
		
		m_nRingXS = m_radiiXS.size();
	}
	
	// Read in the azimuthal subdivisions
	{
		stringstream inBuf(input.child("sub_azi").child_value());
		int azi;
		inBuf >> azi;
		if(inBuf.fail()){
			Error("Improper input to azimuthal subdivisions!");
		}
		m_subAzi.push_back(azi);
		
		// For now im only supporting one entry (same azi for all rings).
		if(m_subAzi.size() > 1){
			Error("Only supporting on azi type for now.");
		}
	}

	// Read in the radial subdivisions
	{
		stringstream inBuf(input.child("sub_radii").child_value());
		while (!inBuf.eof()) {
			int sub;
			inBuf >> sub;
			m_subRad.push_back(sub);
			
			if(inBuf.fail()){
				stringstream msg;
				msg << "Ran into a problem reading radial subdivisions for pin ID="
				    << m_id; 
				Error(msg.str().c_str());
			}
		}
		if(!inBuf.eof()){
			stringstream msg;
			msg << "Dangling data in radial subdivisions for pin ID="
			    << m_id;
		}
		
		// Make sure we have the same number of radial subdivs as rings
		if(m_subRad.size() != m_nRingXS){
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
	for(int ixs=0; ixs<m_nRingXS; ixs++){
		double vn = (m_radiiXS[ixs]*m_radiiXS[ixs] - rxsi*rxsi) / 
		            m_subRad[ixs];
		
		for(int ir=0; ir<m_subRad[ixs]; ir++){
			double r = sqrt(vn + ri*ri);
			m_radii.push_back(r);
			cout << r << endl;
			ri = r;
		}
		rxsi = m_radiiXS[ixs];
	}
	
	return;
	
}