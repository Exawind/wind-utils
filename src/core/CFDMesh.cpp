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

#include "CFDMesh.h"
#include "PerfUtils.h"

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <Ionit_Initializer.h>

#include <limits>
#include <algorithm>

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
{}

CFDMesh::CFDMesh
(
    stk::ParallelMachine& comm,
    const int ndim
) : comm_(comm),
    meta_(ndim),
    bulk_(meta_, comm),
    stkio_(comm)
{}

void CFDMesh::init(stk::io::DatabasePurpose db_purpose)
{
    if (input_db_ == "")
        throw std::runtime_error("CFDMesh::init called on empty database file");

    const std::string timerName("CFDMesh::init_metadata");
    auto timeMon = get_stopwatch(timerName);
    stkio_.add_mesh_database(input_db_, db_purpose);
    stkio_.set_bulk_data(bulk_);
    stkio_.create_input_mesh();
    stkio_.add_all_mesh_fields_as_input_fields();

    // Everyone needs coordinates
    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(coords, meta_.universal_part(), meta_.spatial_dimension(), nullptr);
}

/**
 * \param output_db Filename for the output Exodus database
 * \param time (Optional) time to write (default = 0.0)
 */
void CFDMesh::write_database
(
    std::string output_db,
    double time
)
{
    const std::string timerName("CFDMesh::write_database");
    auto timeMon = get_stopwatch(timerName);
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

/**
 */
void CFDMesh::write_database_with_fields(std::string output_db)
{
    size_t fh = stkio_.create_output_mesh(output_db, stk::io::WRITE_RESTART);
    auto numSteps = stkio_.get_num_time_steps();
    auto times = stkio_.get_time_steps();

    for (int i=0; i<numSteps; i++) {
        stkio_.read_defined_input_fields(i);
        stkio_.begin_output_step(fh, times[i]);
        stkio_.write_defined_output_fields(fh);
        stkio_.end_output_step(fh);
    }

}

BoxType CFDMesh::calc_bounding_box(const stk::mesh::Selector selector, bool verbose)
{
    const std::string timerName("CFDMesh::calc_bounding_box");
    auto timeMon = get_stopwatch(timerName);

    auto ndim = meta_.spatial_dimension();
    std::vector<double> bBoxMin(3, std::numeric_limits<double>::max());
    std::vector<double> bBoxMax(3, std::numeric_limits<double>::lowest());

    if (ndim == 2) {
        bBoxMin[2] = 0.0;
        bBoxMax[2] = 0.0;
    }

    auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, selector);
    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    for (auto b: bkts) {
        double* pt = stk::mesh::field_data(*coords, *b);
        size_t num_nodes = b->size();

        for (size_t in=0; in < num_nodes; in++) {
            for (unsigned int i=0; i<ndim; i++) {
                bBoxMin[i] = std::min(bBoxMin[i], pt[in*ndim + i]);
                bBoxMax[i] = std::max(bBoxMax[i], pt[in*ndim + i]);
            }
        }
    }

    std::vector<double> gMin(3), gMax(3);
    stk::all_reduce_min(
        bulk_.parallel(), bBoxMin.data(), gMin.data(), ndim);
    stk::all_reduce_max(
        bulk_.parallel(), bBoxMax.data(), gMax.data(), ndim);

    PointType minPt(gMin[0], gMin[1], gMin[2]);
    PointType maxPt(gMax[0], gMax[1], gMax[2]);
    BoxType bbox(minPt, maxPt);

    if (verbose && bulk_.parallel_rank() == 0) {
        std::cout << "\nMesh bounding box: \n\t" << bbox;
    }

    return bbox;
}

}  // nalu
}  // sierra
