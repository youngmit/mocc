#pragma once

#include "pugixml.hpp"

#include "source.hpp"

namespace mocc {
	/**
	This returns a unique pointer to the base \ref Source class, which has been
	allocated to the appropriate type for handling the source specified in the
	passed XML tag.

	Right now, only P0 scattering is supported, so a P0 source will always be
	constructed. Using a factory in this context is in anticipation of Pn
	scattering, in which there will be more polymorphism on the source type.
	Even then, it might be desired to treat some aspects of the source as
	template parameters on the sweeper type, or its sweep kernel method,
	somewhat lessening the value of a factory. But we will figure that out when
	we get there.
	*/
    UP_Source_t SourceFactory( const pugi::xml_node &input, 
        int n_reg, const XSMesh *xs_mesh, const ArrayB2 &flux );
}