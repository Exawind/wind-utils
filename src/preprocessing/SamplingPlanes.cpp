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

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>

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
    fluidPart_(),
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
    fluidPart_ = zp["fluid_part"].as<std::string>();
    heights_ = zp["heights"].as<std::vector<double>>();
    name_format_ = zp["part_name_format"].as<std::string>();
    dx_ = zp["dx"].as<double>();
    dy_ = zp["dy"].as<double>();
}

void SamplingPlanes::initialize()
{
    stk::mesh::Part* part = meta_.get_part(fluidPart_);
    if (NULL == part) {
        throw std::runtime_error(
            "SamplingPlanes: Fluid realm not found in mesh database.");
    }

    std::cerr << "SamplingPlanes: Registering parts to meta data:" << std::endl;
    for(auto zh: heights_) {
        std::string pName = (boost::format(name_format_)%zh).str();
        stk::mesh::Part* part_check = meta_.get_part(pName);
        if (part_check != NULL){
            throw std::runtime_error(
                "SamplingPlanes: Cannot overwrite existing part in database: "+ pName);
        } else {
            stk::mesh::Part& part = meta_.declare_part(
                pName, stk::topology::ELEMENT_RANK);
            stk::io::put_io_part_attribute(part);
            stk::mesh::set_topology(part, stk::topology::SHELL_QUAD_4);
            std::cerr << "\t " << pName << std::endl;

            VectorFieldType* coords = meta_.get_field<VectorFieldType>(
                stk::topology::NODE_RANK, "coordinates");
            stk::mesh::put_field(*coords, part, meta_.spatial_dimension());
        }
    }
}

void SamplingPlanes::run()
{
    calc_bounding_box();

    for (auto zh: heights_) {
        generate_zplane(zh);
    }
}

void SamplingPlanes::calc_bounding_box()
{
    stk::mesh::Part* part = meta_.get_part(fluidPart_);
    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    stk::mesh::Selector s_part(*part);
    const stk::mesh::BucketVector& node_buckets = bulk_.get_buckets(
        stk::topology::NODE_RANK, s_part);

    // Setup bounding box array
    for (int i=0; i<ndim_; i++) {
        bBox_[0][i] = std::numeric_limits<double>::max();
        bBox_[1][i] = -std::numeric_limits<double>::max();
    }

    for(size_t ib=0; ib < node_buckets.size(); ib++) {
        stk::mesh::Bucket& bukt = *node_buckets[ib];
        double* pt = stk::mesh::field_data(*coords, bukt);

        for(size_t in=0; in<bukt.size(); in++) {
            for(int i=0; i<ndim_; i++) {
                if (pt[i] < bBox_[0][i]) bBox_[0][i] = pt[i];
                if (pt[i] > bBox_[1][i]) bBox_[1][i] = pt[i];
            }
        }
    }
    std::cerr << "Mesh bounding box: " << std::endl;
    for(size_t i=0; i<2; i++) {
        for(int j=0; j<ndim_; j++)
            std::cerr << "\t" << bBox_[i][j];
        std::cerr << std::endl;
    }

    mx_ = (bBox_[1][0] - bBox_[0][0]) / dx_;
    my_ = (bBox_[1][1] - bBox_[0][1]) / dy_;
    nx_ = mx_ + 1;
    ny_ = my_ + 1;
    std::cerr << "Number of nodes per plane: "
              << (nx_ * ny_) << " [ " << nx_ << " x " << ny_ << " ]"
              << std::endl;
}

void SamplingPlanes::generate_zplane(const double zh)
{
    const std::string pName = (boost::format(name_format_)%zh).str();
    stk::mesh::Part& part = *meta_.get_part(pName);

    unsigned numPoints = nx_ * ny_;
    int numElems = mx_ * my_;
    std::vector<stk::mesh::EntityId> newIDs(numPoints), elemIDs(numElems);
    std::vector<stk::mesh::Entity> nodeVec(numPoints);

    bulk_.modification_begin();
    bulk_.generate_new_ids(stk::topology::NODE_RANK, numPoints, newIDs);
    std::cout << "H = " << zh << "; Part = " << pName
              << "\n\t Node range = [" << newIDs[0]
              << " - " << newIDs[numPoints-1] << "]" << std::endl;
    for(unsigned i=0; i<numPoints; i++) {
        stk::mesh::Entity node = bulk_.declare_entity(
            stk::topology::NODE_RANK, newIDs[i], part);
        nodeVec[i] = node;
    }

    bulk_.generate_new_ids(stk::topology::ELEMENT_RANK, numElems, elemIDs);
    std::cout << "\t Elem range = [" << elemIDs[0] << " - "
              << elemIDs[numElems-1] << "]" << std::endl;

    for(unsigned j=0; j<my_; j++)
        for(unsigned i=0; i<mx_; i++) {
            stk::mesh::EntityIdVector nids(4);
            nids[0] = newIDs[j*nx_ + i];
            nids[1] = newIDs[j*nx_ + i + 1];
            nids[2] = newIDs[(j+1)*nx_ + i + 1];
            nids[3] = newIDs[(j+1)*nx_ + i];
            stk::mesh::declare_element(bulk_, part, elemIDs[j*mx_+i], nids);
        }
    bulk_.modification_end();

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    double xmin = bBox_[0][0];
    double ymin = bBox_[0][1];
    for(unsigned j=0; j<ny_; j++)
        for (unsigned i=0; i<nx_; i++) {
            double* pt = stk::mesh::field_data(*coords, nodeVec[j*nx_+i]);
            pt[0] = xmin + i * dx_;
            pt[1] = ymin + j * dy_;
            pt[2] = zh;
        }
}

}  // nalu
}  // sierra
