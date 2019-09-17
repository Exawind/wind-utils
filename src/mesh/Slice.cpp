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

#include <cmath>

#include "mesh/Slice.h"
#include "core/ParallelInfo.h"
#include "core/PerfUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"

namespace sierra {
namespace nalu {

namespace {
template <typename ViewType>
inline void
cross_prod_helper(const ViewType& v1, const ViewType& v2, ViewType& v3)
{
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

}

Slice::Slice(
    CFDMesh& mesh
) : mesh_(mesh)
{}

void Slice::load(const YAML::Node& node)
{
    int axesLoaded[NDim] = {0, 0, 0};
    // Parse axes
    int numAxes = 0;
    for (int i=0; i < NDim; ++i) {
        std::string axesName = "axis" + std::to_string(i+1);
        if (node[axesName]) {
            axesLoaded[i] = 1;
            const auto axis = node[axesName].as<std::vector<double>>();
            NGP_ThrowRequireMsg(
                axis.size() == NDim, "Incorrect dimensions for axes");

            // Normalize and save the axis vector
            double aMag = 0.0;
            for (int d=0; d < NDim; ++d) {
                axes_(i, d) = axis[d];
                aMag += axis[d] * axis[d];
            }
            aMag = std::sqrt(aMag);

            for (int d=0; d < NDim; ++d)
                axes_(i, d) /= aMag;

            numAxes++;
        }
    }
    NGP_ThrowRequireMsg(
        numAxes > 1, "Need at least two axes provided in the input file");

    for (int i=0; i < NDim; ++i) {
        if (axesLoaded[i] == 1) continue;

        auto tmpvec = Kokkos::subview(axes_, i, Kokkos::ALL);
        cross_prod_helper(
            Kokkos::subview(axes_, (i + 1) % NDim, Kokkos::ALL),
            Kokkos::subview(axes_, (i + 2) % NDim, Kokkos::ALL),
            tmpvec);
    }

    // Load the origin and grid dimensions/information
    const auto origin = node["origin"].as<std::vector<double>>();
    const auto grid_len = node["grid_lengths"].as<std::vector<double>>();
    grid_dx_ = node["grid_dx"].as<std::vector<double>>();
    NGP_ThrowRequireMsg(origin.size() == NDim, "Incorrect dimensions for origin");
    NGP_ThrowRequireMsg(grid_len.size() == 2, "Incorrect dimensions for grid_len");
    NGP_ThrowRequireMsg(grid_dx_.size() == 2, "Incorrect dimensions for grid_dx");

    // Compute the bounding box vertices
    for (int d=0; d < NDim; ++d) {
        vertices_(0, d) = origin[d];
        vertices_(1, d) = origin[d] + axes_(0, d) * grid_len[0];
        vertices_(3, d) = origin[d] + axes_(1, d) * grid_len[1];
        vertices_(2, d) = vertices_(1, d) + vertices_(3, d) - origin[d];
    }

    grid_dims_.resize(2);
    grid_dims_[0] = static_cast<size_t>(grid_len[0] / grid_dx_[0]) + 1;
    grid_dims_[1] = static_cast<size_t>(grid_len[1] / grid_dx_[1]) + 1;

    wind_utils::get_optional(node, "num_planes", num_planes_);
    if (num_planes_ > 1) {
        plane_offsets_ = node["plane_offsets"].as<std::vector<double>>();
        NGP_ThrowRequireMsg(plane_offsets_.size() == num_planes_,
                            "Invalid number of offsets specified");
    } else {
        plane_offsets_.push_back(0.0);
    }

    // Part names that will be generated
    partNamePrefix_ = node["part_name_prefix"].as<std::string>();
}

void Slice::initialize()
{
    const auto& pinfo = get_mpi();

    const std::string timerName("Slice::initialize");
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    auto& coords = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(
        coords, meta.universal_part(), NDim, nullptr);

    pinfo.info()
        << "Slice: Registering parts to meta data: " << std::endl;

    for (int i=0; i < num_planes_; ++i) {
        std::string pName = partNamePrefix_ + "_" + std::to_string(i+1);
        {
            auto* part = meta.get_part(pName);
            NGP_ThrowRequireMsg(
                part == nullptr, "Slice: cannot overwrite an existing part");
        }

        pinfo.info() << "  -  " << pName << std::endl;
        auto& part = meta.declare_part(pName, stk::topology::ELEM_RANK);
        stk::mesh::set_topology(part, stk::topology::SHELL_QUAD_4);
        stk::io::put_io_part_attribute(part);
        stk::mesh::put_field_on_mesh(coords, part, NDim, nullptr);
        partVec_.push_back(&part);
    }
}

void Slice::run()
{
    const auto& pinfo = get_mpi();
    const std::string timerName("Slice::run");
    auto timeMon = get_stopwatch(timerName);

    pinfo.info() << "Generating slices for: " << partNamePrefix_ << std::endl;

    auto& bulk = mesh_.bulk();
    size_t npPlane = grid_dims_[0] * grid_dims_[1];
    size_t numPoints = npPlane * num_planes_;
    size_t num_elems_pp = (grid_dims_[0] - 1) * (grid_dims_[1] - 1);
    size_t num_elems = num_elems_pp * num_planes_;
    std::vector<stk::mesh::EntityId> nodeIDs(numPoints);
    std::vector<stk::mesh::EntityId> elemIDs(num_elems);

    bulk.modification_begin();
    {
        {
            pinfo.info() << "Creating nodes... ";
            bulk.generate_new_ids(stk::topology::NODE_RANK, numPoints, nodeIDs);

            size_t nidx=0;
            unsigned marker = 1;
            for (int i=0; i < num_planes_; ++i) {
                auto* part = partVec_[i];
                for (size_t ip = 0; ip < npPlane; ++ip) {
                    bulk.declare_entity(
                        stk::topology::NODE_RANK, nodeIDs[nidx++], *part);

                    if (marker <= (nidx * 10 / numPoints)) {
                        pinfo.info() << marker * 10 << "% ";
                        marker++;
                    }
                }
            }
            pinfo.info() << std::endl;
        }

        {
            pinfo.info() << "Creating elements... ";
            unsigned marker = 1;
            size_t eidx = 0;
            bulk.generate_new_ids(stk::topology::ELEM_RANK, num_elems, elemIDs);

            stk::mesh::EntityIdVector nids(4);
            for (int i=0; i < num_planes_; ++i) {
                auto* part = partVec_[i];
                size_t offset = i * npPlane;
                for (size_t iy=0; iy < grid_dims_[1]-1; ++iy) {
                    size_t ij = offset + iy * grid_dims_[0];
                    size_t ijp1 = offset + (iy + 1) * grid_dims_[0];
                    for (size_t ix=0; ix < grid_dims_[0]-1; ++ix) {
                        nids[0] = nodeIDs[ij + ix];
                        nids[1] = nodeIDs[ij + ix + 1];
                        nids[2] = nodeIDs[ijp1 + ix + 1];
                        nids[3] = nodeIDs[ijp1 + ix];

                        stk::mesh::declare_element(
                            bulk, *part, elemIDs[eidx++], nids);

                        if (marker <= (eidx * 10) / num_elems) {
                            pinfo.info() << marker * 10 << "% ";
                            marker++;
                        }
                    }
                }
            }
            pinfo.info() << std::endl;
        }
    }
    bulk.modification_end();

    auto* coords = mesh_.meta().get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    pinfo.info() << "Generating coordinate field" << std::endl;
    const double dy = 1.0 / static_cast<double>(grid_dims_[1] - 1);
    const double dx = 1.0 / static_cast<double>(grid_dims_[0] - 1);
    for (int k=0; k < num_planes_; ++k) {
        pinfo.info() << " - " << partVec_[k]->name() << std::endl;
        size_t offset = k * npPlane;
        for (int j=0; j < grid_dims_[1]; ++j) {
            size_t offset1 = offset + j * grid_dims_[0];
            const double ry = j * dy;
            for (int i=0; i < grid_dims_[0]; ++i) {
                size_t idx = offset1 + i;
                const double rx = i * dx;
                auto node = bulk.get_entity(stk::topology::NODE_RANK, nodeIDs[idx]);
                double *pt = stk::mesh::field_data(*coords, node);

                pt[0] = ((1.0 - rx) * (1.0 - ry)) * vertices_(0, 0) +
                    (rx * (1.0 - ry)) * vertices_(1, 0) +
                    (rx * ry) * vertices_(2, 0) +
                    ((1.0 - rx) * ry) * vertices_(3, 0) +
                    plane_offsets_[k] * axes_(2, 0);

                pt[1] = ((1.0 - rx) * (1.0 - ry)) * vertices_(0, 1) +
                    (rx * (1.0 - ry)) * vertices_(1, 1) +
                    (rx * ry) * vertices_(2, 1) +
                    ((1.0 - rx) * ry) * vertices_(3, 1) +
                    plane_offsets_[k] * axes_(2, 1);

                pt[2] = ((1.0 - rx) * (1.0 - ry)) * vertices_(0, 2) +
                    (rx * (1.0 - ry)) * vertices_(1, 2) +
                    (rx * ry) * vertices_(2, 2) +
                    ((1.0 - rx) * ry) * vertices_(3, 2) +
                    plane_offsets_[k] * axes_(2, 2);
            }
        }
    }
}

}  // nalu
}  // sierra
