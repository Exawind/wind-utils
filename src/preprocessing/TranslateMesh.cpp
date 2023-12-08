//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#include "TranslateMesh.h"
#include "core/ClassRegistry.h"
#include "core/PerfUtils.h"

#include "Kokkos_Core.hpp"

#include <vector>
#include <string>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, TranslateMesh, "move_mesh");

TranslateMesh::TranslateMesh(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void
TranslateMesh::load(const YAML::Node& node)
{
    if (node["mesh_parts"]) {
        const auto& fParts = node["mesh_parts"];
        if (fParts.Type() == YAML::NodeType::Scalar) {
            partNames_.push_back(fParts.as<std::string>());
        } else {
            partNames_ = fParts.as<std::vector<std::string>>();
        }
    }

    transVec_ = node["offset_vector"].as<std::vector<double>>();
}

void TranslateMesh::initialize()
{
    auto& meta = mesh_.meta();
    auto nparts = partNames_.size();

    if (nparts > 0) {
        for (auto pName: partNames_) {
            stk::mesh::Part* part = meta.get_part(pName);
            if (NULL == part) {
                throw std::runtime_error(
                    "TranslateMesh:: Part not found in database = " + pName);
            } else {
                parts_.push_back(part);
            }
        }
    } else {
        parts_.push_back(&meta.universal_part());
    }
}

void TranslateMesh::run()
{
    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();

    if (bulk.parallel_rank() == 0)
        std::cout << "Translating mesh" << std::endl;

    const std::string timerName = "TranslateMesh::run()";
    auto timeMon = get_stopwatch(timerName);

    const int ndim = meta.spatial_dimension();
    VectorFieldType* coords = meta.get_field<double>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::Selector s_part = stk::mesh::selectUnion(parts_);
    const auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, s_part);

    using DynamicScheduleType = Kokkos::Schedule<Kokkos::Dynamic>;
    using TeamPolicyType = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, DynamicScheduleType>;
    using TeamHandleType = TeamPolicyType::member_type;

    Kokkos::parallel_for(
        TeamPolicyType(bkts.size(), Kokkos::AUTO),
        [&](const TeamHandleType& team) {
            auto* b = bkts[team.league_rank()];
            double* xyz = stk::mesh::field_data(*coords, *b);

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team, b->size()),
                [&](const size_t& in) {
                    const size_t offset = in * ndim;
                    for (int d = 0; d < ndim; d++)
                        xyz[offset+d] += transVec_[d];
                });
        });

    mesh_.set_write_flag();
}

}  // nalu
}  // sierra
