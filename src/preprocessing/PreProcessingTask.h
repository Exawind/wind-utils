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
