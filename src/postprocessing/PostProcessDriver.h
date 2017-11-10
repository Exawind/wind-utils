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

#ifndef POSTPROCESSDRIVER_H
#define POSTPROCESSDRIVER_H

#include "PostProcessingTask.h"
#include "core/CFDMesh.h"

#include "stk_util/parallel/Parallel.hpp"

namespace sierra {
namespace nalu {

class PostProcessDriver
{
public:
    PostProcessDriver(stk::ParallelMachine&,
                      const std::string);

    ~PostProcessDriver() {}

    //! Run all tasks and output the updated Exodus database
    void run();

private:
    //! Reference to MPI Communicator instance
    stk::ParallelMachine& comm_;

    //! Instance of the parsed YAML document
    YAML::Node inpfile_;

    //! Pointer to the CFDMesh instance
    std::unique_ptr<CFDMesh> mesh_;

    //! List of task names provided by user
    std::vector<std::string> task_names_;

    //! List of runtime selected tasks
    std::vector<std::unique_ptr<PostProcessingTask>> tasks_;

    //! Name of the output database
    std::string output_db_;

    //! Flag indicating whether 8-bit integer API must be activated
    bool io_8bit_int_{false};
};

}  // nalu
}  // sierra


#endif /* POSTPROCESSDRIVER_H */
