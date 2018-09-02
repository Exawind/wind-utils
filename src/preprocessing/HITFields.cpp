//  Copyright 2016 National Renewhite Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applichite law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#include "HITFields.h"
#include "core/YamlUtils.h"
#include "core/KokkosWrappers.h"
#include "core/PerfUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include <fstream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, HITFields, "init_hit_fields");

HITFields::HITFields(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void HITFields::load(const YAML::Node& node)
{
    // Setup mean velocity field
    auto mvel = node["mean_velocity"].as<std::vector<double>>();
    if (mvel.size() !=3)
        throw std::runtime_error("Invalid mean velocity field provided");
    mean_vel_ = mvel;

    // Process part info
    auto fluid_partnames = node["fluid_parts"].as<std::vector<std::string>>();
    fluid_parts_.resize(fluid_partnames.size());

    auto& meta = mesh_.meta();
    for(size_t i=0; i < fluid_partnames.size(); i++) {
        auto* part = meta.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }

    // Get the HIT filename
    hit_filename_ = node["hit_file"].as<std::string>();
    // Get the mesh dimensions
    hit_mesh_dims_ = node["hit_dims"].as<std::vector<int>>();
}

void HITFields::initialize()
{
    const std::string timerName = "HITields::initialize";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    VectorFieldType& velocity = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");
    for (auto part: fluid_parts_) {
        stk::mesh::put_field_on_mesh(velocity, *part, nullptr);
    }
    mesh_.add_output_field("velocity");
}


void HITFields::run()
{
    const std::string timerName = "HITields::run";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();
    const int nDim = meta.spatial_dimension();
    const int nx = hit_mesh_dims_[0];
    const int ny = hit_mesh_dims_[1];
    const int nz = hit_mesh_dims_[2];

    VectorFieldType* velocity = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    std::ifstream hitfile(hit_filename_, std::ios::in | std::ios::binary);
    if (!hitfile.is_open())
        throw std::runtime_error("HITFields:: Error opening file: " + hit_filename_);

    size_t numNodes = (nx * ny * nz);
    size_t numBytes = numNodes * sizeof(double) * 6;
    std::vector<double> buffer(numBytes);

    stk::mesh::Selector sel = stk::mesh::selectUnion(fluid_parts_);
    auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);
    hitfile.read(reinterpret_cast<char*>(buffer.data()), numBytes);

    std::vector<double> minVel = {1.0e10, 1.0e10, 1.0e10};
    std::vector<double> maxVel = {-1.0e10, -1.0e10, -1.0e10};
    for (size_t ib=0; ib < bkts.size(); ib++) {
      auto& b = *bkts[ib];

      for (size_t in=0; in < b.size(); in++) {
        auto node = b[in];
        auto nodeID = bulk.identifier(node);
        // Determine the offset into the buffer array
        //
        // Assume the same ordering of mesh nodes as in the HIT file
        //
        // Skip the (x, y, z) entries for this nodes
        size_t idx = get_index(nodeID) * 6 + 3;
        double* vel = stk::mesh::field_data(*velocity, node);

        for (int d=0; d < nDim; d++) {
          vel[d] = mean_vel_[d] + buffer[idx + d];
          minVel[d] = std::min(vel[d], minVel[d]);
          maxVel[d] = std::max(vel[d], maxVel[d]);
        }
      }
    }
    for (int d=0; d < nDim; d++)
      std::cout << "    Vel[" << d << "]: min = "
                << minVel[d] << "; max = " << maxVel[d] << std::endl;
}

size_t HITFields::get_index(size_t nodeid)
{
    const size_t nx = hit_mesh_dims_[0];
    const size_t ny = hit_mesh_dims_[1];
    const size_t nz = hit_mesh_dims_[2];

    const size_t nxny = (nx + 1) * (ny + 1);
    const size_t nx1 = (nx + 1);

    size_t nid0 = nodeid - 1;
    const size_t iz = (nid0 / nxny) % nz;
    nid0 %= nxny;

    const size_t iy = (nid0 / nx1) % ny;
    const size_t ix = nid0 % nx1 % nx;

    return (iz * (nx * ny) + iy * nx + ix);
}

}  // nalu
}  // sierra
