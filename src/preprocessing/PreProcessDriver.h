#ifndef PREPROCESSDRIVER_H
#define PREPROCESSDRIVER_H

#include "PreProcessingTask.h"
#include "core/CFDMesh.h"

#include "stk_util/parallel/Parallel.hpp"

#include <memory>
#include <vector>
#include <string>

namespace sierra {
namespace nalu {

class PreProcessDriver
{
public:
    PreProcessDriver(stk::ParallelMachine&,
                     const std::string);

    ~PreProcessDriver() {}

    void run();

private:
    stk::ParallelMachine& comm_;

    YAML::Node inpfile_;

    std::unique_ptr<CFDMesh> mesh_;

    std::vector<std::unique_ptr<PreProcessingTask>> tasks_;

    std::string output_db_;
};

} // nalu
} // sierra

#endif /* PREPROCESSDRIVER_H */
