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

#include "PreProcessDriver.h"
#include "core/YamlUtils.h"
#include <iostream>

namespace sierra {
namespace nalu {

/**
 * \param comm MPI Communicator reference
 * \param filename Name of the YAML input file
 */
PreProcessDriver::PreProcessDriver
(
    stk::ParallelMachine& comm,
    const std::string filename
) : comm_(comm),
    inpfile_(YAML::LoadFile(filename))
{
    const YAML::Node& data = inpfile_["nalu_preprocess"];
    std::string input_db = data["input_db"].as<std::string>();
    output_db_ = data["output_db"].as<std::string>();
    wind_utils::get_optional(data, "transfer_fields", transfer_fields_);

    if (data["ioss_8bit_ints"]) {
        io_8bit_int_ = data["ioss_8bit_ints"].as<bool>();
    }

    // Create a new StkIO instance
    mesh_.reset(new CFDMesh(comm, input_db));

    // Mesh decomposition option
    if (data["automatic_decomposition_type"]) {
        std::string decompType = data["automatic_decomposition_type"].as<std::string>();
        mesh_->stkio().property_add(
            Ioss::Property("DECOMPOSITION_METHOD", decompType));
    }
    // 8-bit integer fixes
    if (io_8bit_int_) {
        mesh_->stkio().property_add(Ioss::Property("INTEGER_SIZE_DB",8));
        mesh_->stkio().property_add(Ioss::Property("INTEGER_SIZE_API",8));
    }
    // Initialize the Mesh MetaData
    mesh_->init();

    task_names_ = data["tasks"].as<std::vector<std::string>>();

    tasks_.resize(task_names_.size());
    for (size_t i=0; i < task_names_.size(); i++) {
        tasks_[i].reset(PreProcessingTask::create(*mesh_, data, task_names_[i]));
    }

    if (stk::parallel_machine_rank(comm) == 0) {
        std::cout << "Found " << task_names_.size() << " tasks\n";
        for (auto ts: task_names_) {
            std::cout << "    - " << ts << "\n";
        }
        std::cout << std::endl;
    }
}

void PreProcessDriver::run()
{
    bool dowrite = (stk::parallel_machine_rank(comm_) == 0);

    // Perform metadata updates
    if (dowrite) std::cout << "Performing metadata updates... " << std::endl;
    for (auto& t: tasks_)
        t->initialize();
    if (dowrite) std::cout << "Metadata update completed" << std::endl;

    // Load the bulk data
    if (dowrite) std::cout << "Reading mesh bulk data... ";
    mesh_->stkio().populate_bulk_data();
    if (dowrite) std::cout << "done." << std::endl;

    // Perform modifications to bulk data
    for (size_t i=0; i < tasks_.size(); i++) {
        auto& t = tasks_[i];
        if (dowrite)
            std::cout
                << "\n--------------------------------------------------\n"
                << "Begin task: " << task_names_[i] << std::endl;
        t->run();
        if (dowrite) std::cout << "End task: " << task_names_[i] << std::endl;
     }

    stk::parallel_machine_barrier(mesh_->bulk().parallel());
    if (stk::parallel_machine_rank(comm_) == 0)
        std::cout << "\nAll tasks completed; writing mesh... " << std::endl;

    if (transfer_fields_)
        mesh_->write_database_with_fields(output_db_);
    else
        mesh_->write_database(output_db_);
    if (stk::parallel_machine_rank(comm_) == 0)
        std::cout << "Exodus results file: " << output_db_ << std::endl;
}

} // nalu
} // sierra
