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


#include "ABLFields.h"
#include "core/LinearInterpolation.h"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, ABLFields, "init_abl_fields");

ABLFields::ABLFields(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    fluid_parts_(0),
    ndim_(meta_.spatial_dimension()),
    doVelocity_(false),
    doTemperature_(false)
{
    load(node);
}

void ABLFields::load(const YAML::Node& abl)
{
    auto fluid_partnames = abl["fluid_parts"].as<std::vector<std::string>>();

    if (abl["velocity"]) {
        doVelocity_ = true;
        load_velocity_info(abl["velocity"]);
    }

    if (abl["temperature"]) {
        doTemperature_ = true;
        load_temperature_info(abl["temperature"]);
    }

    fluid_parts_.resize(fluid_partnames.size());

    for(size_t i=0; i<fluid_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }
}

void ABLFields::initialize()
{
    if (doVelocity_) {
        VectorFieldType& velocity = meta_.declare_field<VectorFieldType>(
            stk::topology::NODE_RANK, "velocity");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field(velocity, *part);
        }
        mesh_.add_output_field("velocity");
    }

    if (doTemperature_) {
        ScalarFieldType& temperature = meta_.declare_field<ScalarFieldType>(
            stk::topology::NODE_RANK, "temperature");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field(temperature, *part);
        }
        mesh_.add_output_field("temperature");
    }
}

void ABLFields::run()
{
    if (bulk_.parallel_rank() == 0)
        std::cerr << "Generating ABL fields" << std::endl;
    if (doVelocity_) init_velocity_field();

    if (doTemperature_) init_temperature_field();
}

void ABLFields::load_velocity_info(const YAML::Node& abl)
{
    vHeights_ = abl["heights"].as<std::vector<double>>();
    auto nHeights = vHeights_.size();

    auto velInputs = abl["values"].as<std::vector<std::vector<double>>>();
    ThrowAssertMsg(
        (nHeights == velInputs.size()),
        "ABLFields: Mismatch between sizes of heights and velocities provided "
        "for initializing ABL fields. Check input file.");

    ThrowAssertMsg(
        (ndim_ == velInputs.at(0).size()),
        "ABLFields: Velocity components have all 3 components");

    velocity_.resize(ndim_);
    for (int i=0; i<ndim_; i++) {
        velocity_[i].resize(nHeights);
    }

    // Transpose the arrays for ease of interpolation of components
    for (size_t i=0; i<nHeights; i++)
        for (int j=0; j<ndim_; j++) {
            velocity_[j][i] = velInputs[i][j];
        }
}

void ABLFields::load_temperature_info(const YAML::Node& abl)
{
    THeights_ = abl["heights"].as<std::vector<double>>();
    TValues_ = abl["values"].as<std::vector<double>>();

    ThrowAssertMsg(
        (THeights_.size() == TValues_.size()),
        "ABLFields: Mismatch between sizes of heights and temperature values provided"
        "for initializing ABL fields. Check input file.");
}

void ABLFields::init_velocity_field()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* vel = stk::mesh::field_data(*velocity, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
            const double zh = xyz[in*ndim_ + 2];

            for (int j=0; j<ndim_; j++) {
                utils::linear_interp(
                    vHeights_, velocity_[j], zh, vel[in * ndim_ + j]);
            }
        }
    }
}

void ABLFields::init_temperature_field()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* temperature = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "temperature");

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* temp = stk::mesh::field_data(*temperature, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
            const double zh = xyz[in*ndim_ + 2];
            utils::linear_interp(THeights_, TValues_, zh, temp[in]);
        }
    }
}


} // nalu
} // sierra
