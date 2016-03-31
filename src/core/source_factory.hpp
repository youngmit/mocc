/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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