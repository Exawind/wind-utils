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
#include "core/YamlUtils.h"
#include "core/KokkosWrappers.h"
#include "core/PerfUtils.h"

#include <random>

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
    periodicParts_(0),
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
    const std::string timerName = "ABLfields::initialize";
    auto timeMon = get_stopwatch(timerName);
    if (doVelocity_) {
        VectorFieldType& velocity = meta_.declare_field<VectorFieldType>(
            stk::topology::NODE_RANK, "velocity");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field_on_mesh(velocity, *part, nullptr);
        }
        mesh_.add_output_field("velocity");
    }

    if (doTemperature_) {
        ScalarFieldType& temperature = meta_.declare_field<ScalarFieldType>(
            stk::topology::NODE_RANK, "temperature");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field_on_mesh(temperature, *part, nullptr);
        }
        mesh_.add_output_field("temperature");
    }
}

void ABLFields::run()
{
    const std::string timerName = "ABLfields::run";
    auto timeMon = get_stopwatch(timerName);
    if (bulk_.parallel_rank() == 0)
        std::cout << "Generating ABL fields" << std::endl;
    if (doVelocity_) init_velocity_field();

    if (doTemperature_) init_temperature_field();

    mesh_.set_write_flag();
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

    if (abl["perturbations"]) {
        perturbU_ = true;
        auto& pnode = abl["perturbations"];
        zRefHeight_ = pnode["reference_height"].as<double>();
        std::vector<double> amplitude(2, 0.0);
        amplitude = pnode["amplitude"].as<std::vector<double>>();
        if (amplitude.size() != 2)
            throw std::runtime_error("ABLFields: Invalid size for velocity perturbation amplitude array.");
        deltaU_ = amplitude[0];
        deltaV_ = amplitude[1];
        std::vector<double> periods(2, 0.0);
        periods = pnode["periods"].as<std::vector<double>>();
        if (periods.size() != 2)
            throw std::runtime_error("ABLFields: Invalid size for velocity perturbation periods array.");
        Uperiods_ = periods[0];
        Vperiods_ = periods[1];
    }
}

void ABLFields::load_temperature_info(const YAML::Node& abl)
{
    THeights_ = abl["heights"].as<std::vector<double>>();
    TValues_ = abl["values"].as<std::vector<double>>();

    if (abl["perturbations"]) {
        perturbT_ = true;
        auto& pnode = abl["perturbations"];
        thetaCutoffHt_ = pnode["cutoff_height"].as<double>();
        thetaAmplitude_ = pnode["amplitude"].as<double>();

        if (pnode["skip_periodic_parts"])
            periodicParts_ = pnode["skip_periodic_parts"].as<std::vector<std::string>>();

        wind_utils::get_optional(pnode, "random_gauss_mean", thetaGaussMean_);
        wind_utils::get_optional(pnode, "random_gauss_var", thetaGaussVar_);
    }

    ThrowAssertMsg(
        (THeights_.size() == TValues_.size()),
        "ABLFields: Mismatch between sizes of heights and temperature values provided"
        "for initializing ABL fields. Check input file.");
}

void ABLFields::init_velocity_field()
{
    const std::string timerName = "ABLfields::init_velocity_field";
    auto timeMon = get_stopwatch(timerName);
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    Kokkos::parallel_for(
        TeamPolicyType(fluid_bkts.size(), Kokkos::AUTO),
        [&](const TeamMemberType &team) {
            stk::mesh::Bucket &fbkt = *fluid_bkts[team.league_rank()];
            double *xyz = stk::mesh::field_data(*coords, fbkt);
            double *vel = stk::mesh::field_data(*velocity, fbkt);

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team, fbkt.size()),
                [&](const size_t& in) {
                    const double zh = xyz[in * ndim_ + 2];

                    for (int j = 0; j < ndim_; j++) {
                        utils::linear_interp(vHeights_, velocity_[j], zh,
                                             vel[in * ndim_ + j]);
                    }
                });
        });

    if (perturbU_) perturb_velocity_field();
}

void ABLFields::perturb_velocity_field()
{
    /** The velocity perturbations are similar to the streaks that are added to
     * the velocity field in NREL's SOWFA solution.
     *
     */
    // TODO: Handle height as a distance from terrain

    const std::string timerName = "ABLfields::perturb_velocity_field";
    auto timeMon = get_stopwatch(timerName);
    const double pi = std::acos(-1.0);
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    auto bbox = mesh_.calc_bounding_box(fluid_union, false);
    const double xMin = bbox.get_x_min();
    const double xMax = bbox.get_x_max();
    const double yMin = bbox.get_y_min();
    const double yMax = bbox.get_y_max();

    const double aval = (Uperiods_ * 2.0 * pi / (yMax - yMin));
    const double bval = (Vperiods_ * 2.0 * pi / (xMax - xMin));
    const double ufac = deltaU_ * std::exp(0.5) / zRefHeight_;
    const double vfac = deltaV_ * std::exp(0.5) / zRefHeight_;

    Kokkos::parallel_for(
        TeamPolicyType(fluid_bkts.size(), Kokkos::AUTO),
        [&](const TeamMemberType &team) {
            stk::mesh::Bucket& fbkt = *fluid_bkts[team.league_rank()];
            double* xyz = stk::mesh::field_data(*coords, fbkt);
            double* vel = stk::mesh::field_data(*velocity, fbkt);

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team, fbkt.size()),
                [&](const size_t& in) {
                    const size_t offset = in * ndim_;
                    const double zh = xyz[offset + 2];

                    const double xl = xyz[offset] - xMin;
                    const double yl = xyz[offset + 1] - yMin;
                    const double zl = (zh / zRefHeight_);
                    const double damp = std::exp(-0.5 * zl * zl);

                    // U perturbations
                    vel[offset] += ufac * damp * zh * std::cos(aval * yl);
                    // V perturbations
                    vel[offset + 1] += vfac * damp * zh * std::sin(bval * xl);
                    // No perturbations for w component
                });
        });
}

void ABLFields::init_temperature_field()
{
    const std::string timerName = "ABLfields::init_temperature_field";
    auto timeMon = get_stopwatch(timerName);
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* temperature = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "temperature");

    Kokkos::parallel_for(
        TeamPolicyType(fluid_bkts.size(), Kokkos::AUTO),
        [&](const TeamMemberType &team) {
            stk::mesh::Bucket& fbkt = *fluid_bkts[team.league_rank()];
            double* xyz = stk::mesh::field_data(*coords, fbkt);
            double* temp = stk::mesh::field_data(*temperature, fbkt);

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(team, fbkt.size()),
                [&](const size_t& in) {

                    const double zh = xyz[in*ndim_ + 2];
                    utils::linear_interp(THeights_, TValues_, zh, temp[in]);
                });
        });

    if (perturbT_) perturb_temperature_field();
}

void ABLFields::perturb_temperature_field()
{
    /** Perturbations for the temperature field is based on the following paper:
     *
     *  D. Munoz-Esparza, B. Kosovic, J. van Beeck, J. D. Mirocha, A stocastic
     *  perturbation method to generate inflow turbulence in large-eddy
     *  simulation models: Application to neutrally stratified atmospheric
     *  boundary layers. Physics of Fluids, Vol. 27, 2015.
     *
     */

    const std::string timerName = "ABLfields::perturb_temperature_field";
    auto timeMon = get_stopwatch(timerName);
    stk::mesh::PartVector partVec(0);
    for (size_t i=0; i < periodicParts_.size(); i++) {
        auto* part = meta_.get_part(periodicParts_[i]);
        if (part != nullptr)
            partVec.push_back(part);
    }

    stk::mesh::Selector fluid_union = (
        stk::mesh::selectUnion(fluid_parts_)
        & !stk::mesh::selectUnion(partVec));
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* temperature = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "temperature");

    // Random number generator for adding temperature perturbations
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> gaussNormal{thetaGaussMean_, thetaGaussVar_};

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* temp = stk::mesh::field_data(*temperature, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
            const double zh = xyz[in*ndim_ + 2];

            if (zh < thetaCutoffHt_)
                temp[in] += thetaAmplitude_ * gaussNormal(gen);
        }
    }
}

} // nalu
} // sierra
