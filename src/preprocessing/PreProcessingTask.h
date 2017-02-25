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

#ifndef PREPROCESSINGTASK_H
#define PREPROCESSINGTASK_H

#include "core/ClassRegistry.h"
#include "core/CFDMesh.h"

#include "yaml-cpp/yaml.h"

namespace sierra {
namespace nalu {

/** An abstract implementation of a PreProcessingTask
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
class PreProcessingTask
{
public:
    /**
     * \param mesh A CFDMesh instance
     */
    PreProcessingTask(CFDMesh& mesh) : mesh_(mesh) {}

    virtual ~PreProcessingTask() {}

    virtual void initialize() = 0;

    virtual void run() = 0;

    DECLARE_INHERITANCE_REGISTRY
    (
        PreProcessingTask,
        (
            CFDMesh& mesh,
            const YAML::Node& node
        ),
        (mesh, node)
    );

    /** Runtime creation of concrete task instance
     */
    static PreProcessingTask* create(
        CFDMesh&,
        const YAML::Node&,
        std::string);

protected:
    //! Reference to the CFDMesh instance
    CFDMesh& mesh_;

private:
    PreProcessingTask() = delete;
    PreProcessingTask(const PreProcessingTask&) = delete;
};

}  // nalu
}  // sierra

#endif /* PREPROCESSINGTASK_H */
