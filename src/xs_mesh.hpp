#pragma once

#include <vector>

#include "core_mesh.hpp"

namespace mocc {
    class XSMesh {
    public:
        // Default constructor does almost nothing, and lets some other code
        // tell it what to do
        XSMesh() {}

        // XSMesh provides its own facility to initialize itself from a CoreMesh
        XSMesh( const CoreMesh& mesh );

        // Return the number of energy groups
        unsigned int n_grp() const {
            return ng_;
        }
        

    private:
        // list of mesh regions corresponding to each XSMesh region
        std::vector<VecI> reg_;
        VecF xsmactr_;
        VecF xsmacnf_;
        VecF xsmackf_;
        VecF xsmacch_;
        unsigned int ng_;
    };
}
