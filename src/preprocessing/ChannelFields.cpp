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
    doTKE_(false),
    length_(0.0),
    width_(0.0),
    height_(0.0),
    C_(0.0),
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

    if (channel["turbulent_ke"].as<bool>()) {
        doTKE_ = true;
        load_tke_info(channel["tke"]);
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
            stk::mesh::put_field_on_mesh(velocity, *part, nullptr);
        }
        mesh_.add_output_field("velocity");
    }

    if (doTKE_) {
        ScalarFieldType& tke = meta_.declare_field<ScalarFieldType>(
            stk::topology::NODE_RANK, "turbulent_ke");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field_on_mesh(tke, *part, nullptr);
        }
        mesh_.add_output_field("turbulent_ke");
    }
}

void ChannelFields::run()
{
    if (bulk_.parallel_rank() == 0)
        std::cout << "Generating channel fields" << std::endl;
    setup_parameters();
    if (doVelocity_) init_velocity_field();
    if (doTKE_) init_tke_field();

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

void ChannelFields::load_tke_info(const YAML::Node&)
{}

void ChannelFields::setup_parameters()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    auto bbox = mesh_.calc_bounding_box(fluid_union,false);
    length_ = bbox.get_x_max() - bbox.get_x_min();
    height_ = bbox.get_y_max() - bbox.get_y_min();
    width_ = bbox.get_z_max() - bbox.get_z_min();
    const double delta = 0.5 * height_;
    utau_ = Re_tau_ * viscosity_ / delta;
    const double yph = height_ * utau_ / viscosity_;
    C_ = (1.0/kappa_ * log(1.0 + kappa_ * yph)) / (1 - exp(- yph / 11.0) - yph / 11.0 * exp(- yph / 3)) + log(kappa_) / kappa_;
}

double ChannelFields::reichardt(const double y)
{
    const double yp = std::min(y, height_ - y) * utau_ / viscosity_;
    return (1.0/kappa_ * log(1.0 + kappa_ * yp)) + (C_ - log(kappa_) / kappa_) * (1 - exp(- yp / 11.0) - yp / 11.0 * exp(- yp / 3));
}

double ChannelFields::u_perturbation(const double x, const double, const double z)
{
    const double pert_u = a_pert_u_ * sin(k_pert_u_ * M_PI / width_ * z) * sin(k_pert_u_ * M_PI / length_ * x);
    const double rand_u = a_rand_u_ * (2. * (double)rand() / RAND_MAX - 1);
    return pert_u * rand_u;
}

double ChannelFields::w_perturbation(const double x, const double, const double z)
{
    const double pert_w = a_pert_w_ * sin(k_pert_w_ * M_PI / length_ * x)  * sin(k_pert_w_ * M_PI / width_ * z);
    const double rand_w = a_rand_w_ * (2. * (double)rand() / RAND_MAX - 1);
    return pert_w * rand_w;
}

void ChannelFields::init_velocity_field()
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
	  const double x = xyz[in*ndim_ + 0];
	  const double y = xyz[in*ndim_ + 1];
	  const double z = xyz[in*ndim_ + 2];

	  const double pert_u = u_perturbation(x,y,z);
	  const double pert_w = w_perturbation(x,y,z);

	  vel[in * ndim_ + 0] = utau_ * (reichardt(y) + pert_u);
	  vel[in * ndim_ + 1] = 0.0;
	  vel[in * ndim_ + 2] = utau_ * pert_w;
        }
    }
}

void ChannelFields::init_tke_field()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* tke = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "turbulent_ke");

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
    	double* tke_tmp = stk::mesh::field_data(*tke, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
    	  const double y = xyz[in*ndim_ + 1];

	  tke_tmp[in] = 0.5 * utau_ * utau_ * reichardt(y);
        }
    }
}


} // nalu
} // sierra
