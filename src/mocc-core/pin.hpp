#pragma once

#include <vector>
#include <memory>
#include <map>

#include "pin_mesh.hpp"
#include "global_config.hpp"


namespace mocc {
	// The Pin class is a concrete instantiaion of a physical pin. It
	// essentially applies materials to regions of a PinMesh. Nothin' fancy.
	class Pin {
	public:
		Pin( const pugi::xml_node &input, 
            const std::map<int, UP_PinMesh_t> &meshes );
        ~Pin(){
            return;
        }

        const PinMesh& mesh() const {
            return *pin_mesh_;
        }

        int id() const {
            return id_;
        }

        int n_reg() const {
            return pin_mesh_->n_reg();
        }

        int mesh_id() const {
            return pin_mesh_->id();
        }

        real_t vol() const {
            return pin_mesh_->vol();
        }

        const VecF& vols() const {
            return pin_mesh_->vols();
        }

        const VecI& mat_ids() const {
            return mat_IDs_;
        }
	private:
        // Pin ID
        const unsigned int id_;
        unsigned int mesh_id_;
		// Immutable reference to the pin mesh object (owned by CoreMesh)
		PinMesh const * pin_mesh_;
		// Material IDs to apply to each XS region of the pin mesh
		VecI mat_IDs_;
	};

    typedef std::shared_ptr<Pin> SP_Pin_t;
    typedef std::unique_ptr<Pin> UP_Pin_t;
}
