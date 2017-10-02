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


#include "ChannelFields.h"
#include "core/LinearInterpolation.h"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, ChannelFields, "init_channel_fields");

ChannelFields::ChannelFields(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    fluid_parts_(0),
    ndim_(meta_.spatial_dimension()),
    doVelocity_(false),
    Re_tau_(0.0),
    viscosity_(0.0)
{
    load(node);
    srand(seed_);
}

void ChannelFields::load(const YAML::Node& channel)
{
    auto fluid_partnames = channel["fluid_parts"].as<std::vector<std::string>>();

    if (channel["velocity"]) {
        doVelocity_ = true;
        load_velocity_info(channel["velocity"]);
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

void ChannelFields::initialize()
{
    if (doVelocity_) {
        VectorFieldType& velocity = meta_.declare_field<VectorFieldType>(
            stk::topology::NODE_RANK, "velocity");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field(velocity, *part);
        }
        mesh_.add_output_field("velocity");
    }
}

void ChannelFields::run()
{
    if (bulk_.parallel_rank() == 0)
        std::cout << "Generating channel fields" << std::endl;
    if (doVelocity_) init_velocity_field();

    mesh_.set_write_flag();
}

void ChannelFields::load_velocity_info(const YAML::Node& channel)
{
  if (channel["Re_tau"])
    Re_tau_ = channel["Re_tau"].as<double>();
  else
    throw std::runtime_error("ChannelFields: missing mandatory Re_tau parameter");

  if (channel["viscosity"])
    viscosity_ = channel["viscosity"].as<double>();
  else
    throw std::runtime_error("ChannelFields: missing mandatory viscosity parameter");
}

void ChannelFields::init_velocity_field()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    auto bbox = mesh_.calc_bounding_box(fluid_union,false);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    const double length = bbox.get_x_max() - bbox.get_x_min();
    const double height = bbox.get_y_max() - bbox.get_y_min();
    const double width = bbox.get_z_max() - bbox.get_z_min();
    const double delta = 0.5 * height;
    const double utau = Re_tau_ * viscosity_ / delta;
    const double yph = height * utau / viscosity_;
    const double C = (1.0/kappa_ * log(1.0 + kappa_ * yph)) / (1 - exp(- yph / 11.0) - yph / 11.0 * exp(- yph / 3)) + log(kappa_) / kappa_;

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* vel = stk::mesh::field_data(*velocity, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
	  const double x = xyz[in*ndim_ + 0];
	  const double y = xyz[in*ndim_ + 1];
	  const double z = xyz[in*ndim_ + 2];

	  const double yp = std::min(y, height - y) * utau / viscosity_;
	  const double reichardt = (1.0/kappa_ * log(1.0 + kappa_ * yp)) + (C - log(kappa_) / kappa_) * (1 - exp(- yp / 11.0) - yp / 11.0 * exp(- yp / 3));

	  const double pert_u = sin(k_pert_u_ * M_PI / width * z);
	  const double pert_w = sin(k_pert_w_ * M_PI / length * x);

	  const double rand_u = a_rand_u_ * 2. * (double)rand() / RAND_MAX - 1;
	  const double rand_w = a_rand_w_ * 2. * (double)rand() / RAND_MAX - 1;

	  vel[in * ndim_ + 0] = utau * (reichardt + pert_u + rand_u);
	  vel[in * ndim_ + 1] = 0.0;
	  vel[in * ndim_ + 2] = utau * (pert_w + rand_w);
        }
    }
}


} // nalu
} // sierra
