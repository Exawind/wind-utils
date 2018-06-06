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

#ifndef HITFIELDS_H
#define HITFIELDS_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

class HITFields: public PreProcessingTask
{
public:
    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    HITFields(CFDMesh&, const YAML::Node&);

    virtual ~HITFields() {}

    //! Declare velocity and temperature fields and register them for output
    void initialize();

    //! Initialize the velocity and/or temperature fields by linear interpolation
    void run();

private:
    HITFields() = delete;
    HITFields(const HITFields&) = delete;

    size_t get_index(size_t);

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    stk::mesh::PartVector fluid_parts_;

    std::vector<double> mean_vel_{0.0, 0.0, 0.0};

    std::vector<int> hit_mesh_dims_{0, 0, 0};

    std::string hit_filename_;

};

}  // nalu
}  // sierra


#endif /* HITFIELDS_H */
