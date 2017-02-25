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

/** A driver that runs all preprocessor tasks.
 *
 * This class is responsible for reading the input file, parsing the
 * user-requested list of tasks, initializing the task instances, executing
 * them, and finally writing out the updated Exodus database with changed
 * inputs.
 */
class PreProcessDriver
{
public:
    PreProcessDriver(stk::ParallelMachine&,
                     const std::string);

    ~PreProcessDriver() {}

    //! Run all tasks and output the updated Exodus database
    void run();

private:
    //! Reference to MPI Communicator instance
    stk::ParallelMachine& comm_;

    //! Instance of the parsed YAML document
    YAML::Node inpfile_;

    //! Pointer to the CFDMesh instance
    std::unique_ptr<CFDMesh> mesh_;

    //! List of runtime selected tasks
    std::vector<std::unique_ptr<PreProcessingTask>> tasks_;

    //! Name of the output database
    std::string output_db_;
};

} // nalu
} // sierra

#endif /* PREPROCESSDRIVER_H */
