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

#include "ABLStatistics.h"

#include "stk_util/parallel/ParallelReduce.hpp"

#include <cmath>
#include <fstream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PostProcessingTask, ABLStatistics, "abl_statistics");

ABLStatistics::ABLStatistics
(
    CFDMesh& mesh,
    const YAML::Node& node
) : PostProcessingTask(mesh),
    meta_(mesh_.meta()),
    bulk_(mesh_.bulk()),
    fluid_parts_(0),
    ndim_(meta_.spatial_dimension())
{
    if (ndim_ != 3)
        throw std::runtime_error("ABLStatistics is only supported for 3-D meshes");

    // Load some default field names
    field_map_["velocity"] = "velocity";
    field_map_["temperature"] = "temperature";
    field_map_["sfs_stress"] = "sfs_stress";
    field_map_["temperature_resolved_stress"] = "temperature_resolved_stress";
    field_map_["temperature_variance"] = "temperature_variance";
    load(node);
}

void ABLStatistics::load(const YAML::Node& node)
{
    bool dowrite = (bulk_.parallel_rank() == 0);
    auto fluid_partnames = node["fluid_parts"].as<std::vector<std::string>>();

    fluid_parts_.resize(fluid_partnames.size());
    size_t ii = 0;
    for (auto name: fluid_partnames) {
        stk::mesh::Part* part = meta_.get_part(name);
        if (nullptr == part) {
            throw std::runtime_error("ABLStatistics:: Missing fluid part in mesh database: " +
                                     name);
        } else {
            fluid_parts_[ii++] = part;
        }
    }

    auto& field_name_map = node["field_map"];
    for (auto it: field_name_map)
        field_map_[it.first.as<std::string>()] = it.second.as<std::string>();

    // Implement constant spacing only for now
    auto& htnode = node["height_info"];
    zmin_ = htnode["min_height"].as<double>();
    zmax_ = htnode["max_height"].as<double>();
    dz_ = htnode["delta_height"].as<double>();

    nheights_ = static_cast<int>((zmax_ - zmin_) / dz_) + 1;
    heights_.resize(nheights_);
    for (int i=0; i < nheights_; i++)
        heights_[i] = zmin_ + i * dz_;

    node_counters_.resize(nheights_);
    velMean_.resize(nheights_ * ndim_);
    tempMean_.resize(nheights_);
    sfsMean_.resize(nheights_ * ndim_ * 2);

    if (dowrite) {
        std::cerr << "ABLStatistics:: Will process data at " << nheights_ << " heights." << std::endl;
    }
}

void
ABLStatistics::initialize()
{
    auto& velocity = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, field_map_["velocity"]);
    auto& temperature = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, field_map_["temperature"]);
    stk::mesh::FieldBase *sfs_stress = &meta_.declare_field<
        stk::mesh::Field<double, stk::mesh::SimpleArrayTag>>(
        stk::topology::NODE_RANK, field_map_["sfs_stress"]);

    for (auto* part: fluid_parts_) {
        stk::mesh::put_field_on_mesh(velocity, *part, ndim_, nullptr);
        stk::mesh::put_field_on_mesh(temperature, *part, 1, nullptr);
        stk::mesh::put_field_on_mesh(*sfs_stress, *part, ndim_*2, nullptr);
    }
}

void
ABLStatistics::populate_solution()
{
    bool dowrite = (bulk_.parallel_rank() == 0);
    auto num_steps = mesh_.stkio().get_num_time_steps();
    auto times = mesh_.stkio().get_time_steps();

    auto final_time = times[num_steps - 1];
    std::vector<stk::io::MeshField> missing_fields;
    auto found_time = mesh_.stkio().read_defined_input_fields(final_time, &missing_fields);

    if (missing_fields.size() > 0) {
        if (dowrite) {
            std::cout << "Missing fields in the solution file: " << std::endl;
            for (size_t i=0; i < missing_fields.size(); i++)
                std::cout << "    -" << missing_fields[i].field()->name() << std::endl;
        }
        throw std::runtime_error("ABLStatistics:: missing required fields in database");
    }

    if (dowrite)
        std::cout << "ABLStatistics:: Found " << num_steps << " time steps in database. Using time = "
                  << found_time << " to calculate ABL statistics" << std::endl;
}

void
ABLStatistics::average_planes()
{
    const stk::mesh::Selector sel = stk::mesh::selectUnion(fluid_parts_);
    auto bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    const VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    const VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, field_map_["velocity"]);
    const ScalarFieldType* temperature = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, field_map_["temperature"]);

    // Initialize mean arrays to zero before we start accumulation
    for (int ih=0; ih < nheights_; ih++) {
        tempMean_[ih] = 0.0;

        for (int d=0; d < ndim_; d++) {
            velMean_[ih * ndim_ + d] = 0.0;
        }
    }

    // Process node buckets in this MPI rank and sum up velocities on a plane
    for (auto b: bkts) {
        for (size_t in =0; in < b->size(); in++) {
            auto node = (*b)[in];

            // Determine the index from the z-coordinate. Here we are assuming
            // constant z-spacing as well as flat terrain.
            double* crd = stk::mesh::field_data(*coords, node);
            int ih = static_cast<int>(std::floor(((crd[2] - zmin_) + 1.0e-10) / dz_));
            node_counters_[ih] += 1;

            // Velocity calculations
            double* vel_vec = stk::mesh::field_data(*velocity, node);
            for (int d=0; d < ndim_; d++)
                velMean_[ih * ndim_ + d] += vel_vec[d];

            // temperature calculations
            double* temp = stk::mesh::field_data(*temperature, node);
            tempMean_[ih] += *temp;

            // TODO: Add SFS stress calculations here
        }
    }

    // Perform global sum and average
    std::vector<double> gVelMean(nheights_ * ndim_, 0.0);
    std::vector<double> gTempMean(nheights_, 0.0);
    std::vector<int> gNodeCtr(nheights_, 0);

    stk::all_reduce_sum(bulk_.parallel(), velMean_.data(), gVelMean.data(), nheights_ * ndim_);
    stk::all_reduce_sum(bulk_.parallel(), tempMean_.data(), gTempMean.data(), nheights_);
    stk::all_reduce_sum(bulk_.parallel(), node_counters_.data(), gNodeCtr.data(), nheights_);

    // Reassign averages to member data structures
    for (int ih=0; ih < nheights_; ih++) {
        tempMean_[ih] = gTempMean[ih] / gNodeCtr[ih];

        for (int d=0; d < ndim_; d++)  {
            velMean_[ih * ndim_ + d] = gVelMean[ih * ndim_ + d] / gNodeCtr[ih];
        }
    }
}

void
ABLStatistics::output_averages()
{
    // Only output data files on MPI master rank
    if (bulk_.parallel_rank() != 0) return;

    std::ofstream velfile;
    std::ofstream tempfile;

    // TODO: Provide user defined output filenames
    velfile.open("abl_velocity_stats.dat", std::ofstream::out);
    tempfile.open("abl_temperature_stats.dat", std::ofstream::out);

    velfile << "# Height, Ux, Uy, Uz" << std::endl;
    tempfile << "# Height, T" << std::endl;

    // TODO: Provide precision changes
    for (int ih=0; ih < nheights_; ih++) {
        tempfile << heights_[ih] << " " << tempMean_[ih] << std::endl;

        velfile << heights_[ih];
        for (int d=0; d < ndim_; d++)
            velfile << " " << velMean_[ih * ndim_ + d];
        velfile << std::endl;
    }
    velfile.close();
    tempfile.close();

    std::cout << "ABLStatistics:: Finished writing statistics files" << std::endl;
}

void
ABLStatistics::run()
{
    // Load the solution from the database
    populate_solution();

    // Process nodes on mesh and compute spatial averages for the fields
    average_planes();

    // Output to text files
    output_averages();
}

}  // nalu
}  // sierra
