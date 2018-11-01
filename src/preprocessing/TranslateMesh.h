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

#ifndef TRANSLATEMESH_H
#define TRANSLATEMESH_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/** Translate a mesh
 */
class TranslateMesh: public PreProcessingTask
{
public:
    TranslateMesh(CFDMesh&, const YAML::Node&);

    virtual ~TranslateMesh() {}

    virtual void initialize();

    virtual void run();

private:
    TranslateMesh() = delete;
    TranslateMesh(const TranslateMesh&) = delete;

    void load(const YAML::Node&);

    //! Part names of the mesh that need to be translated
    std::vector<std::string> partNames_;

    //! Parts of the mesh that needs to be translated
    stk::mesh::PartVector parts_;

    //! Translation vector
    std::vector<double> transVec_;
};

}  // nalu
}  // sierra


#endif /* TRANSLATEMESH_H */
