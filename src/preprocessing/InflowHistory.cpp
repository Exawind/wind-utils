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

#include "preprocessing/InflowHistory.h"

#include "core/YamlUtils.h"
#include "core/KokkosWrappers.h"
#include "core/PerfUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include <fstream>
#include <vector>
#include <string>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, InflowHistory, "time_varying_inflow");

InflowHistory::InflowHistory(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void InflowHistory::load(const YAML::Node& node)
{
    auto partnames = node["fluid_parts"].as<std::vector<std::string>>();
    partVec_.resize(partnames.size());

    auto& meta = mesh_.meta();
    for(size_t i=0; i < partnames.size(); i++) {
        auto* part = meta.get_part(partnames[i]);
        if (nullptr == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     partnames[i]);
        } else {
            partVec_[i] = part;
        }
    }

    inflow_filename_ = node["inflow_file"].as<std::string>();
    output_db_ = node["time_history_db"].as<std::string>();
    numSteps_ = node["num_timesteps"].as<int>();
}

void InflowHistory::initialize()
{
    const std::string timerName = "InflowHistory::initialize";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    VectorFieldType& velocity = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");
    for (auto part: partVec_) {
        stk::mesh::put_field_on_mesh(
            velocity, *part, meta.spatial_dimension(), nullptr);
    }
    mesh_.add_output_field("velocity");
}

void InflowHistory::run()
{
    const std::string timerName = "InflowHistory::run";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();
    const int nDim = meta.spatial_dimension();

    const stk::mesh::Selector sel = stk::mesh::selectUnion(partVec_);
    auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);

    std::ifstream inflow(inflow_filename_, std::ios::in);
    if (!inflow.is_open())
        throw std::runtime_error(
            "InflowHistory:: Error opening file: " + inflow_filename_);

    auto* velocity = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    double time, uvel, vvel, wvel;
    mesh_.write_timesteps(
        output_db_, numSteps_,
        [&](int tstep) {
            inflow >> time >> uvel >> vvel >> wvel;

            for (auto b: bkts) {
                for (auto node: *b) {
                    double* vel = stk::mesh::field_data(*velocity, node);
                    vel[0] = uvel;
                    vel[1] = vvel;
                    vel[2] = wvel;
                }
            }
            return time;
        });
}

}  // nalu
}  // sierra
