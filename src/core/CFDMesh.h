#ifndef CFDMESH_H
#define CFDMESH_H

#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Field.hpp"
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <string>
#include <memory>
#include <unordered_set>

namespace sierra {
namespace nalu {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;

class CFDMesh
{
public:
    CFDMesh(stk::ParallelMachine&,
            const std::string);

    ~CFDMesh() {}

    inline stk::mesh::MetaData& meta() { return meta_; }

    inline stk::mesh::BulkData& bulk() { return bulk_; }

    inline stk::io::StkMeshIoBroker& stkio() { return stkio_; }

    /** Register a field for output during write
     */
    inline void add_output_field(const std::string field)
    {
        output_fields_.insert(field);
    }

    /** Write the Exodus results database with modifications
     */
    void write_database(std::string output_db, double time=0.0);

private:
    CFDMesh() = delete;
    CFDMesh(const CFDMesh&) = delete;

    //! Instance of the parallel MPI_COMM object
    stk::ParallelMachine& comm_;

    //! STK Metadata object
    stk::mesh::MetaData meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData bulk_;

    //! Filename for the mesh input database
    std::string input_db_;

    //! STK Mesh I/O handler
    stk::io::StkMeshIoBroker stkio_;

    //! List of fields registered for output
    std::unordered_set<std::string> output_fields_;
};

}  // nalu
}  // sierra

#endif /* CFDMESH_H */
