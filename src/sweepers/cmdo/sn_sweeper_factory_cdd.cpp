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

#include "sn_sweeper_factory_cdd.hpp"
#include "util/string_utils.hpp"
#include "sn_sweeper_factory.hpp"

namespace mocc {
namespace cmdo {
CDDPair_t SnSweeperFactory_CDD(const pugi::xml_node &input,
                               const CoreMesh &mesh)
{
    std::string equation = "cdd";
    if (!input.attribute("equation").empty()) {
        equation = input.attribute("equation").value();
        sanitize(equation);
    }

    std::unique_ptr<SnSweeper> sweeper;
    std::shared_ptr<CorrectionData> corrections;

    if (equation != "cdd") {
        LogScreen
            << "Something wants a CDD sweeper, but the equation "
               "specified is different. Keep in mind that the correction data "
               "generated here probably isn't being used"
            << std::endl;
        // \todo this is a cyclical dependency. While not illegal, figure a
        // way around it. Maybe dont actually call the CDD factory from the
        // one called below.
        sweeper     = SnSweeperFactory(input, mesh);
        corrections = std::make_shared<CorrectionData>(
            mesh, sweeper->ang_quad().ndir() / 2, sweeper->n_group());
        corrections->from_data(input);
        return CDDPair_t(std::move(sweeper), corrections);
    }

    // Determine the type of axial treatment and create the right type of
    // sweeper.
    std::string axial = "dd";
    if (!input.attribute("axial").empty()) {
        axial = input.attribute("axial").value();
        sanitize(axial);
    }
    if (axial == "dd") {
        LogScreen << "Diamond Difference axial treatment" << std::endl;
        cmdo::SnSweeper_CDD<CellWorker_CDD_DD> *swp =
            new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_DD>(input, mesh);
        corrections = std::make_shared<CorrectionData>(
            mesh, swp->ang_quad().ndir() / 2, swp->n_group());
        if (!input.child("data").empty()) {
            corrections->from_data(input);
        }

        swp->set_corrections(corrections);
        sweeper.reset(swp);
    } else if (axial == "sc") {
        LogScreen << "Step Characteristics axial treatment" << std::endl;
        cmdo::SnSweeper_CDD<CellWorker_CDD_SC> *swp =
            new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_SC>(input, mesh);
        corrections = std::make_shared<CorrectionData>(
            mesh, swp->ang_quad().ndir() / 2, swp->n_group());
        if (!input.child("data").empty()) {
            corrections->from_data(input);
        }

        swp->set_corrections(corrections);
        sweeper.reset(swp);
    } else if (axial == "fw") {
        LogScreen << "Forward Difference axial treatment" << std::endl;
        cmdo::SnSweeper_CDD<CellWorker_CDD_FW> *swp =
            new cmdo::SnSweeper_CDD<cmdo::CellWorker_CDD_FW>(input, mesh);

        corrections = std::make_shared<CorrectionData>(
            mesh, swp->ang_quad().ndir() / 2, swp->n_group());
        if (!input.child("data").empty()) {
            corrections->from_data(input);
        }

        swp->set_corrections(corrections);
        sweeper.reset(swp);
    } else {
        throw EXCEPT("Unrecognized axial treatment in CDD.");
    }
    return CDDPair_t(std::move(sweeper), corrections);
}
}
}
