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

#ifndef POSTPROCESSINGTASK_H
#define POSTPROCESSINGTASK_H

#include "core/ClassRegistry.h"
#include "core/CFDMesh.h"

#include "core/YamlUtils.h"

namespace sierra {
namespace nalu {

/** An abstract implementation of a PostProcessingTask
 *
 *  This class defines the interface for a pre-processing task and contains the
 *  infrastructure to allow concrete implementations of pre-processing tasks to
 *  register themselves for automatic runtime discovery. Derived classes must
 *  implement two methods:
 *
 * - `initialize` - Perform actions on STK MetaData before processing BulkData
 *
 * - `run` - All actions on BulkData and other operations on mesh after it has
 *    been loaded from the disk.
 *
 * For automatic class registration, the derived classes must implement a
 * constructor that takes two arguments: a CFDMesh reference, and a `const`
 * reference to YAML::Node that contains the inputs necessary for the concrete
 * task implementation. It is the derived class' responsibility to process the
 * input dictionary and perform error checking. No STK mesh manipulations must
 * occur in the constructor.
 *
 */
class PostProcessingTask
{
public:
    PostProcessingTask(CFDMesh& mesh) : mesh_(mesh) {}

    virtual ~PostProcessingTask() {}

    virtual void initialize() = 0;

    virtual void run() = 0;


    DECLARE_INHERITANCE_REGISTRY
    (
        PostProcessingTask,
        (
            CFDMesh& mesh,
            const YAML::Node& node
        ),
        (mesh, node)
    );

    /** Runtime creation of concrete task instance
     */
    static PostProcessingTask* create(
        CFDMesh&,
        const YAML::Node&,
        std::string);

protected:
    //! Reference to the CFDMesh instance
    CFDMesh& mesh_;

private:
    PostProcessingTask() = delete;
    PostProcessingTask(const PostProcessingTask&) = delete;
};

}  // nalu
}  // sierra

#endif /* POSTPROCESSINGTASK_H */
