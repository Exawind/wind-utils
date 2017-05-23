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
    std::string task_type = lookup;
    if (inp["task_type"]) {
        task_type = inp["task_type"].as<std::string>();
    }
    auto it = PreProcessingTaskReg_ConstructorTable_->find(task_type);
    if (it != PreProcessingTaskReg_ConstructorTable_->end()) {
        return (it->second)(mesh, inp);
    } else {
        std::cerr << "ERROR: Invalid PreProcessingTask => " << task_type << std::endl;
        std::cerr << "Valid task types are: " << std::endl;
        for (const auto& t: *PreProcessingTaskReg_ConstructorTable_) {
            std::cerr << "\t" << t.first << std::endl;
        }
    }
    return nullptr;
}

}  // nalu
}  // sierra
