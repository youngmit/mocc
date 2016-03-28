#pragma once

#include <memory>
#include <utility>

#include "pugixml.hpp"

#include "core/core_mesh.hpp"

#include "correction_data.hpp"
#include "sn_sweeper_cdd.hpp"


namespace mocc { namespace cmdo {
    /**
     * Generate an \ref SnSweeper_CDD and associated correction data based on
     * the input provided.
     *
     * This factory is responsible for interpreting the provided input to
     * determine and create the appropriate \ref sn::SnSweeper and return it.
     * This factory is distinct from the vanilla \ref SnSweeperFactory in that
     * it also creates \ref CorrectionData for the sweeper and returns it as
     * well.  This is important, because the sweeper that is ultimately returned
     * from this function is of the base type \ref sn::SnSweeper, which doesn't
     * actually know anything about the existence of correction factors.
     *
     * \note While it would be possible to maintain type information about the
     * CDD nature for the returned sweeper, it would be necessary to propagate
     * the template parameter as well, which becomes unwieldy when a sweeper
     * ends up owning an \ref SnSweeper_CDD as a member. In this case it becomes
     * necessary to template that sweeper class as well.
     *
     * \note This method has one big potential gotcha; the \c std::pair that is
     * returned from the factory, and \ref CorrectionData that it contains is
     * the \a only reference to the \ref CorrectionData that survives the to
     * the end of this function. It is also impossible to get a new one, since
     * sufficient type information to do so is discarded when an \ref
     * sn::SnSweeper is returned. Moral of the story is to be careful with what
     * you do with the return value of this function. Here is a simple example
     * of how to screw up:
     * \code
     * std::unique_ptr<SnSweeper> sweeper = SnSweeperFactory_CDD( input, mesh).first;
     * std::shared_ptr<CorrectionData> corrections = SnSweeperFactory_CDD( input, mesh ).second;
     * \endcode
     *
     */
    CDDPair_t SnSweeperFactory_CDD( const pugi::xml_node &input,
                                    const CoreMesh &mesh);
} }
