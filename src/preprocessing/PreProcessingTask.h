#ifndef PREPROCESSINGTASK_H
#define PREPROCESSINGTASK_H

#include "core/ClassRegistry.h"
#include "core/CFDMesh.h"

#include "yaml-cpp/yaml.h"

namespace sierra {
namespace nalu {

class PreProcessingTask
{
public:
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

    static PreProcessingTask* create(
        CFDMesh&,
        const YAML::Node&,
        std::string);

protected:
    CFDMesh& mesh_;

private:
    PreProcessingTask() = delete;
    PreProcessingTask(const PreProcessingTask&) = delete;
};

}  // nalu
}  // sierra

#endif /* PREPROCESSINGTASK_H */
