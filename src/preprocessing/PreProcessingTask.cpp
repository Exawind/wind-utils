
#include "PreProcessingTask.h"

#include <iostream>

namespace sierra {
namespace nalu {

DEFINE_INHERITANCE_REGISTRY(PreProcessingTask);

PreProcessingTask*
PreProcessingTask::create(
    CFDMesh& mesh,
    const YAML::Node& node,
    std::string lookup)
{
    if (!node[lookup]) {
        throw std::runtime_error("Cannot find input section for task: " + lookup);
    }

    const YAML::Node& inp = node[lookup];
    auto it = PreProcessingTaskReg_ConstructorTable_.find(lookup);
    if (it != PreProcessingTaskReg_ConstructorTable_.end()) {
        return (it->second)(mesh, inp);
    } else {
        std::cerr << "ERROR: Invalid PreProcessingTask => " << lookup << std::endl;
        std::cerr << "Valid task types are: " << std::endl;
        for (const auto& t: PreProcessingTaskReg_ConstructorTable_) {
            std::cerr << "\t" << t.first << std::endl;
        }
    }
    return nullptr;
}

}  // nalu
}  // sierra
