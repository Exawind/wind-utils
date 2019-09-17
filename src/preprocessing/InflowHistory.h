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

#ifndef INFLOWHISTORY_H
#define INFLOWHISTORY_H

#include "preprocessing/PreProcessingTask.h"

namespace sierra {
namespace nalu {

class InflowHistory: public PreProcessingTask
{
public:
    InflowHistory(CFDMesh&, const YAML::Node&);

    virtual ~InflowHistory() = default;

    void initialize();

    void run();

private:
    void load(const YAML::Node&);

    stk::mesh::PartVector partVec_;

    std::string inflow_filename_;

    std::string output_db_{"inflow_history.exo"};

    size_t numSteps_{0};
};

}  // nalu
}  // sierra


#endif /* INFLOWHISTORY_H */
