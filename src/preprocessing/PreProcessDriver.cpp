
#include "PreProcessDriver.h"

namespace sierra {
namespace nalu {

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

    mesh_.reset(new CFDMesh(comm, input_db));

    auto task_names = data["tasks"].as<std::vector<std::string>>();

    tasks_.resize(task_names.size());
    for (size_t i=0; i < task_names.size(); i++) {
        tasks_[i].reset(PreProcessingTask::create(*mesh_, data, task_names[i]));
    }

    std::cerr << "    Found " << task_names.size() << " tasks\n";
    for (auto ts: task_names) {
        std::cerr << "        - " << ts << "\n";
    }
    std::cerr << std::endl;
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

    std::cerr << "\nAll tasks completed; writing mesh... " << std::endl;
    mesh_->write_database(output_db_);
    std::cerr << "Exodus results file: " << output_db_ << std::endl;
}

} // nalu
} // sierra
