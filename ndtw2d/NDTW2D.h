#ifndef NDTW2D_H
#define NDTW2D_H

#include "stk_util/parallel/Parallel.hpp"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Field.hpp"

#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include "yaml-cpp/yaml.h"

#include <vector>

namespace sierra {
namespace nalu {

class NDTW2D
{
public:
    NDTW2D(
        stk::ParallelMachine&,
        const YAML::Node&);

    /**
     * Generate wall distance and output new mesh database
     */
    void run();

private:
    NDTW2D() = delete;
    NDTW2D(const NDTW2D&) = delete;

    //! Load user inputs from YAML input file
    void load(const YAML::Node&);

    //! Generate meta data
    void generate_meta();

    //! Calculate wall distance
    void calc_ndtw();

    //! Spatial dimension (currently only 2)
    const int ndim_;

    //! Instance of parallel MPI_COMM for use with StkMeshIoBroker
    stk::ParallelMachine& comm_;

    //! STK Metadata object
    stk::mesh::MetaData meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData bulk_;

    //! Filename for the input mesh database
    std::string input_db_;

    //! Filename for the output mesh database
    std::string output_db_;

    //! Field name for wall distance 
    std::string wall_dist_name_;

    //! Parts of the fluid mesh where wall distance is calculated
    stk::mesh::PartVector fluid_parts_;

    //! Part names of the wall boundaries
    stk::mesh::PartVector wall_parts_;

    //! STK Mesh I/O handler
    std::unique_ptr<stk::io::StkMeshIoBroker> stkIo_;
};


}  // nalu
}  // sierra

#endif /* NDTW2D_H */
