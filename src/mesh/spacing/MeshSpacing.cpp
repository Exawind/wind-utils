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

#include "MeshSpacing.h"

#include <iostream>

namespace sierra {
namespace nalu {

DEFINE_INHERITANCE_REGISTRY(MeshSpacing);

MeshSpacing*
MeshSpacing::create(
    int npts,
    const YAML::Node& node,
    std::string lookup)
{
    auto it = MeshSpacingReg_ConstructorTable_->find(lookup);
    if (it != MeshSpacingReg_ConstructorTable_->end()) {
        return (it->second)(npts, node);
    } else {
        std::cout << "ERROR: Invalid mesh spacing => " << lookup << std::endl;
        std::cout << "Valid spacing types are: " << std::endl;
        for (const auto& t: *MeshSpacingReg_ConstructorTable_) {
            std::cout << "\t" << t.first << std::endl;
        }
        throw std::runtime_error("Invalid mesh spacing specified.");
    }

    return nullptr;
}

}  // nalu
}  // sierra
