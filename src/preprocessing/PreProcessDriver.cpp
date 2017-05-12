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

    auto task_names = data["tasks"].as<std::vector<std::string>>();

    tasks_.resize(task_names.size());
    for (size_t i=0; i < task_names.size(); i++) {
        tasks_[i].reset(PreProcessingTask::create(*mesh_, data, task_names[i]));
    }

    if (stk::parallel_machine_rank(comm) == 0) {
        std::cerr << "    Found " << task_names.size() << " tasks\n";
        for (auto ts: task_names) {
            std::cerr << "        - " << ts << "\n";
        }
        std::cerr << std::endl;
    }
}

void PreProcessDriver::run()
{
    // Perform metadata updates
    for (auto& t: tasks_)
        t->initialize();

    // Load the bulk data
    mesh_->stkio().populate_bulk_data();

    // Perform modifications to bulk data
    for (auto& t: tasks_)
        t->run();

    stk::parallel_machine_barrier(mesh_->bulk().parallel());
    if (stk::parallel_machine_rank(comm_) == 0)
        std::cerr << "\nAll tasks completed; writing mesh... " << std::endl;
    mesh_->write_database(output_db_);
    if (stk::parallel_machine_rank(comm_) == 0)
        std::cerr << "Exodus results file: " << output_db_ << std::endl;
}

} // nalu
} // sierra
