
#include "CFDMesh.h"

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>

namespace sierra {
namespace nalu {

CFDMesh::CFDMesh
(
    stk::ParallelMachine& comm,
    const std::string filename
) : comm_(comm),
    meta_(),
    bulk_(meta_, comm),
    input_db_(filename),
    stkio_(comm)
{
    stkio_.add_mesh_database(input_db_, stk::io::READ_MESH);
    stkio_.set_bulk_data(bulk_);
    stkio_.create_input_mesh();
    stkio_.add_all_mesh_fields_as_input_fields();

    // Everyone needs coordinates
    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(coords, meta_.universal_part(), meta_.spatial_dimension());
}

void CFDMesh::write_database
(
    std::string output_db,
    double time
)
{
    size_t fh = stkio_.create_output_mesh(output_db, stk::io::WRITE_RESULTS);
    for (auto fname: output_fields_ ) {
        stk::mesh::FieldBase* fld = stk::mesh::get_field_by_name(fname, meta_);
        if (fld != NULL) {
            stkio_.add_field(fh, *fld, fname);
        }
    }
    stkio_.begin_output_step(fh, time);
    stkio_.write_defined_output_fields(fh);
    stkio_.end_output_step(fh);
}

}  // nalu
}  // sierra
