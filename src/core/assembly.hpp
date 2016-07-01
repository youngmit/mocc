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

#include <cassert>
#include <map>
#include <memory>
#include <vector>
#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "lattice.hpp"

namespace mocc {
class Assembly {
public:
    Assembly(const pugi::xml_node &input,
             const std::map<int, UP_Lattice_t> &lattices);

    ~Assembly();

    int id() const
    {
        return id_;
    }

    /*
     * Return the number of pins along the x dimension
     */
    unsigned int nx() const
    {
        return lattices_[0]->nx();
    }

    /**
     * Return the number of pins along the y dimension
     */
    unsigned int ny() const
    {
        return lattices_[0]->ny();
    }

    /**
     * Return the number of planes
     */
    unsigned int nz() const
    {
        return nz_;
    }

    /**
     * Return the total height of the indexed plane in the \ref Assembly
     */
    real_t dz(unsigned int iz) const
    {
        return dz_[iz];
    }

    /**
     * Return a const reference to the vector of plane heights.
     */
    const VecF &dz() const
    {
        return dz_;
    }

    // Return the total size of the assembly in the x dimension
    real_t hx() const
    {
        return hx_;
    }

    // Return the total size of the assembly in the y dimension
    real_t hy() const
    {
        return hy_;
    }

    // Return the total number of FSRs in the assembly
    unsigned int n_reg() const
    {
        return n_reg_;
    }

    // Return the total number of XS regions in the assembly
    unsigned int n_xsreg() const
    {
        return n_xsreg_;
    }

    /**
     * Return a reference to the indexed \ref Lattice.
     */
    const Lattice &operator[](unsigned int iz) const
    {
        assert((iz >= 0) & (iz < lattices_.size()));
        return *lattices_[iz];
    }

    /**
     * \brief Return whether the passed \ref Assembly is compatible with this
     * \ref Assembly
     *
     * In this context, "compatible" means that the assemblies are the same
     * height, have the same plane heights throughout, and have the same
     * subplane parameters throughout.
     */
    bool compatible(const Assembly &other) const;

    /**
     * \brief Return the subplane parameters
     *
     * Subplane parameters are a sequence of integers, each representing the
     * number of planes to be bound into a "macroplane," from bottom to top.
     */
    const VecI &subplane() const
    {
        return subplane_;
    }

private:
    int id_;
    unsigned int nz_;
    VecF dz_;

    real_t hx_;
    real_t hy_;

    size_t n_reg_;
    size_t n_xsreg_;

    // Sub-plane factors. List of numbers of planes that should be bound
    // together, from bottom to top
    VecI subplane_;

    std::vector<const Lattice *> lattices_;
};

typedef std::unique_ptr<Assembly> UP_Assembly_t;

std::map<int, UP_Assembly_t>
ParseAssemblies(const pugi::xml_node &input,
                const std::map<int, UP_Lattice_t> lattices);
}
