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

#include "PostProcessDriver.h"
#include "core/YamlUtils.h"

namespace sierra {
namespace nalu {

PostProcessDriver::PostProcessDriver
(
    stk::ParallelMachine& comm,
    const std::string filename
) : comm_(comm),
    inpfile_(YAML::LoadFile(filename))
{
    const auto& data = inpfile_["nalu_postprocess"];
    std::string input_db = data["input_db"].as<std::string>();

    wind_utils::get_optional(data, "ioss_8bit_ints", io_8bit_int_);

    mesh_.reset(new CFDMesh(comm, input_db));

    // Initialize the Mesh MetaData
    mesh_->init(stk::io::READ_RESTART);

    task_names_ = data["tasks"].as<std::vector<std::string>>();

    tasks_.resize(task_names_.size());
    for (size_t i=0; i < task_names_.size(); i++) {
        tasks_[i].reset(PostProcessingTask::create(*mesh_, data, task_names_[i]));
    }

    if (stk::parallel_machine_rank(comm) == 0) {
        std::cout << "Found " << task_names_.size() << " tasks\n";
        for (auto ts: task_names_) {
            std::cout << "    - " << ts << "\n";
        }
        std::cout << std::endl;
    }
}

void
PostProcessDriver::run()
{
    bool dowrite = (stk::parallel_machine_rank(comm_) == 0);

    if (dowrite) std::cout << "Performing metadata updates.. " << std::endl;
    for (auto& t: tasks_) t->initialize();
    if (dowrite) std::cout << "Metadata update completed" << std::endl;

    if (dowrite) std::cout << "Loading mesh bulk data... ";
    mesh_->stkio().populate_bulk_data();
    if (dowrite) std::cout << "done" << std::endl;

    for (size_t i=0; i < tasks_.size(); i++) {
        if (dowrite)
            std::cout << "Begin task: " << task_names_[i] << std::endl;
        tasks_[i]->run();
        if (dowrite)
            std::cout << "End task: " << task_names_[i] << std::endl;
    }

    stk::parallel_machine_barrier(mesh_->bulk().parallel());
    if (mesh_->db_modified()) {
        if (stk::parallel_machine_rank(comm_) == 0)
            std::cout << "\nAll tasks completed; writing mesh... " << std::endl;

        mesh_->write_database(output_db_);
        if (stk::parallel_machine_rank(comm_) == 0)
            std::cout << "Exodus results file: " << output_db_ << std::endl;
    } else {
        if (stk::parallel_machine_rank(comm_) == 0)
            std::cout << "Input mesh DB not modified; skipping write" << std::endl;
    }
}

}  // nalu
}  // sierra
