#pragma once

#include "pin_mesh.hpp"
#include "global_config.hpp"

#include <vector>

namespace mocc {
	// The Pin class is a concrete instantiaion of a physical pin. It
	// essentially applies materials to regions of a PinMesh. Nothin' fancy.
	class Pin {
	public:
		Pin();
	private:
		// Immutable reference to the pin mesh object (owned by CoreMesh)
		PinMesh const * const pin_mesh_;
		// Material IDs to apply to each XS region of the pin mesh
		VecI mat_IDs_;
	};
}