#pragma once

#include "pin_mesh.hpp"
#include "global_config.hpp"

#include <vector>
#include <memory>

namespace mocc {
	// The Pin class is a concrete instantiaion of a physical pin. It
	// essentially applies materials to regions of a PinMesh. Nothin' fancy.
	class Pin {
	public:
		Pin( int id, PinMesh* pin, VecI mat );
        ~Pin(){
            return;
        }

        PinMesh const * const mesh() const {
            return pin_mesh_;
        }
	private:
        // Pin ID
        const int id_;
		// Immutable reference to the pin mesh object (owned by CoreMesh)
		PinMesh const * const pin_mesh_;
		// Material IDs to apply to each XS region of the pin mesh
		const VecI mat_IDs_;
	};

    typedef std::shared_ptr<Pin> SP_Pin_t;
    typedef std::unique_ptr<Pin> UP_Pin_t;
}
