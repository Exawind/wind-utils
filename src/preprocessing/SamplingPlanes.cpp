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

#include "SamplingPlanes.h"
#include "core/ClassRegistry.h"
#include "core/PerfUtils.h"

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>

#include <boost/format.hpp>

#include <iostream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, SamplingPlanes, "generate_planes");

SamplingPlanes::SamplingPlanes(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    heights_(0),
    bBox_(),
    name_format_(),
    fluidPartNames_(),
    fluidParts_(),
    vertices_(4, std::vector<double>(2,0.0)),
    dx_(0.0),
    dy_(0.0),
    nx_(0),
    ny_(0),
    mx_(0),
    my_(0),
    ndim_(meta_.spatial_dimension())
{
    if (ndim_ != 3)
        throw std::runtime_error("SamplingPlanes only available for 3-D meshes");
    load(node);
}

void SamplingPlanes::load(const YAML::Node& zp)
{
    const auto& fParts = zp["fluid_part"];
    if (fParts.Type() == YAML::NodeType::Scalar) {
        fluidPartNames_.push_back(fParts.as<std::string>());
    } else {
        fluidPartNames_ = fParts.as<std::vector<std::string>>();
    }

    if (zp["boundary_type"]) {
        std::string bdy_type_name = zp["boundary_type"].as<std::string>();

        if (bdy_type_name == "bounding_box")
            bdyType_ = BOUND_BOX;
        else if (bdy_type_name == "quad_vertices")
            bdyType_ = QUAD_VERTICES;
        else
            throw std::runtime_error("Bad option specified for boundary type: " +
                                     bdy_type_name);
    } else {
        bdyType_ = BOUND_BOX;
    }

    heights_ = zp["heights"].as<std::vector<double>>();
    name_format_ = zp["part_name_format"].as<std::string>();

    if (bdyType_ == BOUND_BOX) {
        dx_ = zp["dx"].as<double>();
        dy_ = zp["dy"].as<double>();
    } else {
        mx_ = zp["nx"].as<int>();
        my_ = zp["ny"].as<int>();
        nx_ = mx_ + 1;
        ny_ = my_ + 1;

        vertices_ = zp["vertices"].as<std::vector<std::vector<double>>>();
        if (vertices_.size() != 4)
            throw std::runtime_error("Incorrect number of vertices provided. Expected 4.");
        for (size_t i=0; i<vertices_.size(); i++) {
            size_t nl = vertices_[i].size();
            if ((nl < 2) || (nl > 3))
                throw std::runtime_error(
                    "Inconsistent vertices provided. Check input file.");
        }
    }
}

void SamplingPlanes::initialize()
{
    const std::string timerName = "SamplingPlanes::initialize";
    auto timeMon = get_stopwatch(timerName);
    const auto iproc = bulk_.parallel_rank();
    for (auto pName: fluidPartNames_){
        stk::mesh::Part* part = meta_.get_part(pName);
        if (NULL == part) {
            throw std::runtime_error(
                "SamplingPlanes: Fluid realm not found in mesh database.");
        } else {
            fluidParts_.push_back(part);
        }
    }

    if (iproc == 0)
        std::cout << "SamplingPlanes: Registering parts to meta data:" << std::endl;
    for(auto zh: heights_) {
        std::string pName = (boost::format(name_format_)%zh).str();
        stk::mesh::Part* part_check = meta_.get_part(pName);
        if (part_check != NULL){
            throw std::runtime_error(
                "SamplingPlanes: Cannot overwrite existing part in database: "+ pName);
        } else {
            stk::mesh::Part& part = meta_.declare_part(
                pName, stk::topology::NODE_RANK);
            stk::io::put_io_part_attribute(part);
            // stk::mesh::set_topology(part, stk::topology::SHELL_QUAD_4);

            if (iproc == 0) std::cout << "\t " << pName << std::endl;

            VectorFieldType* coords = meta_.get_field<VectorFieldType>(
                stk::topology::NODE_RANK, "coordinates");
            stk::mesh::put_field_on_mesh(*coords, part, meta_.spatial_dimension(), nullptr);
        }
    }
}

void SamplingPlanes::run()
{
    const std::string timerName = "SamplingPlanes::run";
    auto timeMon = get_stopwatch(timerName);
    calc_bounding_box();

    stk::parallel_machine_barrier(bulk_.parallel());
    for (auto zh: heights_) {
        generate_zplane(zh);
    }

    mesh_.set_write_flag();
}

void SamplingPlanes::calc_bounding_box()
{
    const std::string timerName = "SamplingPlanes::calc_bounding_box";
    auto timeMon = get_stopwatch(timerName);
    auto iproc = bulk_.parallel_rank();
    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    stk::mesh::Selector s_part = stk::mesh::selectUnion(fluidParts_);
    const stk::mesh::BucketVector& node_buckets = bulk_.get_buckets(
        stk::topology::NODE_RANK, s_part);

    // Setup bounding box array
    std::vector<double> bBoxMin(ndim_);
    std::vector<double> bBoxMax(ndim_);
    for (int i=0; i<ndim_; i++) {
        bBoxMin[i] = std::numeric_limits<double>::max();
        bBoxMax[i] = -std::numeric_limits<double>::max();
    }

    for(size_t ib=0; ib < node_buckets.size(); ib++) {
        stk::mesh::Bucket& bukt = *node_buckets[ib];
        double* pt = stk::mesh::field_data(*coords, bukt);

        for(size_t in=0; in<bukt.size(); in++) {
            for(int i=0; i<ndim_; i++) {
                if (pt[i] < bBoxMin[i]) bBoxMin[i] = pt[in * ndim_ + i];
                if (pt[i] > bBoxMax[i]) bBoxMax[i] = pt[in * ndim_ + i];
            }
        }
    }
    stk::all_reduce_min(
        bulk_.parallel(), bBoxMin.data(), bBox_[0].data(), ndim_);
    stk::all_reduce_max(
        bulk_.parallel(), bBoxMax.data(), bBox_[1].data(), ndim_);
    if (iproc == 0) {
        std::cout << "Mesh bounding box: " << std::endl;
        for(size_t i=0; i<2; i++) {
            for(int j=0; j<ndim_; j++)
                std::cout << "\t" << bBox_[i][j];
            std::cout << std::endl;
        }
    }

    if (bdyType_ == BOUND_BOX) {
        mx_ = (bBox_[1][0] - bBox_[0][0]) / dx_;
        my_ = (bBox_[1][1] - bBox_[0][1]) / dy_;
        nx_ = mx_ + 1;
        ny_ = my_ + 1;

        // S-W corner (0 vertex)
        vertices_[0][0] = bBox_[0][0];
        vertices_[0][1] = bBox_[0][1];
        // S-E corner (1 vertex)
        vertices_[1][0] = bBox_[1][0];
        vertices_[1][1] = bBox_[0][1];
        // N-E corner (2 vertex)
        vertices_[2][0] = bBox_[1][0];
        vertices_[2][1] = bBox_[1][1];
        // N-W corner (1 vertex)
        vertices_[3][0] = bBox_[0][0];
        vertices_[3][1] = bBox_[1][1];
    }

    // Reset dx and dy for computations
    dx_ = 1.0 / static_cast<double>(mx_);
    dy_ = 1.0 / static_cast<double>(my_);
    if (iproc == 0){
        std::cout << "Number of nodes per plane: "
                  << (nx_ * ny_) << " [ " << nx_ << " x " << ny_ << " ]"
                  << std::endl;
    }
}

void SamplingPlanes::generate_zplane(const double zh)
{
    const std::string timerName = "SamplingPlanes::generate_zplane";
    auto timeMon = get_stopwatch(timerName);
    const unsigned iproc = bulk_.parallel_rank();
    const unsigned nproc = bulk_.parallel_size();
    const std::string pName = (boost::format(name_format_)%zh).str();
    stk::mesh::Part& part = *meta_.get_part(pName);

    unsigned gNumPoints = nx_ * ny_;
    unsigned numPoints = 0;
    unsigned offset = 0;
    if ((gNumPoints < nproc) && (iproc < gNumPoints)) {
        numPoints = 1;
    } else {
        numPoints = gNumPoints / nproc;
        offset = iproc * numPoints;
        unsigned rem = gNumPoints % nproc;

        if ((rem > 0) && (iproc < rem)) numPoints++;
        offset += (iproc < rem)? iproc : rem;
    }

    std::vector<stk::mesh::EntityId> newIDs(numPoints);
    std::vector<stk::mesh::Entity> nodeVec(numPoints);

    bulk_.modification_begin();
    if (numPoints > 0) {
        bulk_.generate_new_ids(stk::topology::NODE_RANK, numPoints, newIDs);

        for(unsigned i=0; i<numPoints; i++) {
            stk::mesh::Entity node = bulk_.declare_entity(
                stk::topology::NODE_RANK, newIDs[i], part);
            nodeVec[i] = node;
        }
    }
    bulk_.modification_end();

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    for (unsigned k=0; k < numPoints; k++) {
        int j = (offset+k) / nx_;
        int i = (offset+k) % nx_;
        double* pt = stk::mesh::field_data(*coords, nodeVec[k]);

        const double rx = i * dx_;
        const double ry = j * dy_;

        pt[0] = ((1.0 - rx) * (1.0 - ry) * vertices_[0][0] +
                 rx * (1.0 - ry) * vertices_[1][0] +
                 rx * ry * vertices_[2][0] +
                 (1.0 - rx) * ry * vertices_[3][0]);
        pt[1] = ((1.0 - rx) * (1.0 - ry) * vertices_[0][1] +
                 rx * (1.0 - ry) * vertices_[1][1] +
                 rx * ry * vertices_[2][1] +
                 (1.0 - rx) * ry * vertices_[3][1]);
        pt[2] = zh;
    }
}

}  // nalu
}  // sierra
